#read in initial diets and plot intended predators for each prey - just uses proportion of diet with no scaling for predator abundance or prey abundance

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
inputsPath<-paste(basePath, "\\inputs\\",sep="")
plotPath<-paste(DIR$'Figures',"Mortality\\", sep="")
groupsDF<-read.csv(paste(this_path,"CRAM_groups.csv",sep="")); ng<-dim(groupsDF)[1]
predators<-as.character(groupsDF$Code[groupsDF$IsPredator==1]); npreds<-length(predators)

dietsOriginal<-read.csv(paste(inputsPath,"supporting\\CRAM_Diets.csv",sep=""))
diets<-read.csv(paste(inputsPath,"biol_prm\\pPREY\\CRAM_biol_pPREY_original.csv", sep=""))
preyCodesOrig<-colnames(dietsOriginal)[7:(7+ng-1)]; preyCodeIndex<-match(groupsDF$Code, preyCodesOrig)
##overwrite diets with values from diestOriginal. Diets has values ready for pPREY (proportions of prey available to pred)
# but want proportions of diets here
for(j in 2:dim(diets)[2]){
  thisVar<-colnames(diets)[j]; preyAge<-get_first_number(thisVar, n=1); predAge<-get_first_number(thisVar, n=2); thisCode<-gsub(paste("pPREY", predAge,preyAge,sep="|"), "", thisVar)
  xx<-grep(thisCode, dietsOriginal$Code); 
  if(!is.na(preyAge) & !is.na(predAge)){
    xxx<-dietsOriginal[xx:(xx+3),]
    predAgeText<-c("Juvenile", "Adult")[predAge]; preyAgeText<-c("Juvenile", "Adult")[preyAge]
  } else{
    xxx<-dietsOriginal[xx,]
    predAgeText<-"All"; preyAgeText<-"All"
  }  
  index<-xxx$Predator.age==predAgeText & xxx$Prey.age==preyAgeText
  thisDiets<-xxx[index,]
  newColumn<-thisDiets[7:dim(dietsOriginal)[2]]
  diets[1:ng,j]<-as.double(newColumn[preyCodeIndex])
}

preys<-diets[,1]; npreys<-length(preys)
getAgeColor<-function(x){
  thisCol<-myGrey
  if(is.na(x)){
    thisCol<-myGrey
  } else{
    if(x==2){thisCol="black"}
  }
  return(thisCol)
}
getAgeText<-function(x){
  thisText<-"Juv"
  if(is.na(x)){
    thisText<-""
  }else{
    if(x==2){thisText="Adult"}
  }
  return(thisText)
}

pdf(paste(plotPath,"AllIntendedDiets.pdf", sep="")); par(mfrow=c(3,2), mar=c(6,4,2,1), oma=c(0,0,0,0))
for(p in 1:npreds){
  thisPred<-predators[p]
  temp<-diets[,c(1,grep(thisPred, colnames(diets)))]
  test<-sum((temp[,-1]), na.rm=TRUE)
  if(test>0){
    preyAges<-unlist(lapply(colnames(temp)[-1], get_first_number, n=1))
    predAges<-unlist(lapply(colnames(temp)[-1], get_first_number, n=2))
    if(is.na(predAges)){
      
    }else{
      for(predAge in c(1,2)){
        predAgeText<-c("Juvenile", "Adult")[predAge]
        thisData<-temp[,c(1,grep(paste(thisPred,predAge,sep=""), colnames(temp)))]
        xx<-rowSums(thisData[,-1]); index<-xx>0
        thisDiet<-thisData[index,]
        if(dim(thisDiet)[1]>0){
          ymax<-max(thisDiet[,-1],na.rm=TRUE)*1.2
          par(lend=1, las=1)
          plot(x=seq(1,dim(thisDiet)[1])-0.1, y=thisDiet[,2], type="h", lwd=5, col=myGrey, xlab="", ylab="Proportion of diet", xaxt="n", xlim=c(0,dim(thisDiet)[1]+0.5), ylim=c(0,ymax))
          mtext(paste(thisPred, predAgeText,sep=" "), side=3, adj=0)
          points(x=seq(1,dim(thisDiet)[1])+0.1, y=thisDiet[,3], type="h", lwd=5, col="black")
          par(las=2)
          axis(at=seq(1,dim(thisDiet)[1]), labels=thisDiet[,1], side=1)
          legend(legend=c("Juvenile", "Adult"), col=c(myGrey, "black"), lwd=5, seg.len=1, bty="n",x="topleft")
        }
      }
    }
  }
}
dev.off()


############################################################
## read in model run and plot realised diets (will be approx)
# see PredatorFocusPlot.R
