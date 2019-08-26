create_vadt_pred_prey <- function(outdir, fgfile, biolprm, ncout, startyear, toutinc, diet = TRUE){
  # outdir=outPath
  # startyear=modelStartYear
  # toutinc=daysTimeStep
  # diet = TRUE
  # ncout<-"output"
  # fgfile<-paste(thisPath,"..\\CRAM_Groups.csv",sep="")
  # biolprm<-paste(thisPath,"..\\CRAM_base_biol",this_out,".prm",sep="")

  # contants
  nsecs <- 86400
  ndays <- 365
  g_per_ton <- 1e6
  species <- c("BIRD", "FISH", "MAMMAL", "SHARK")
  tons <- function(mgN) return(mgN * 5.7 * 20 / 1e9)
  
  cat("### ------------ Reading in data                                         ------------ ###\n")
  nc_out <- nc_open(paste(outdir, ncout, ".nc", sep = ""))
  # prod_out <- nc_open(paste(outdir, ncout, "PROD.nc", sep = ""))
  # bio_agg <- read.table(paste(outdir, ncout, "BoxBiomass.txt", sep = ""), header = T)
  # ssb <- read.table(paste(outdir, ncout, "SSB.txt", sep = ""), header = TRUE)
  # yoy <- read.table(paste(outdir, ncout, "YOY.txt", sep = ""), header = TRUE)
  bgm <- readLines(paste(outdir, grep(".bgm",dir(outdir), value = T), sep = ""))
  # biomass <- read.table(paste(outdir, ncout, "BiomIndx.txt", sep = ""), header = T)
  # rel_bio <- biomass[,c(1, grep("Rel",colnames(biomass)))]
  # tot_bio <- biomass[,c(1:(grep("Rel",colnames(biomass))[1]-1))]
  fun_group <- read.csv(fgfile, header = T, stringsAsFactors = FALSE)#[, c(1,3, 4,5,6, 9,16, 12)]
  fun_group <- subset(fun_group, fun_group$IsTurnedOn == 1)
  
  biolprm <- readLines(biolprm)
  
    # diet_l <- NULL
    # tot_pred <- NULL
  
  getAverageSN<-function(code,age="both"){
    spp <- fun_group[fun_group$Code == code,c("Name")]
    spp <- str_trim(spp)
    if(age=="Juvenile"){
      SNnames<-paste(spp,c(" Juvenile"),sep="")
    }else if(age=="Adult"){
      SNnames<-paste(spp,c(" Adult"),sep="")
    }else{
      SNnames<-paste(spp,c(" Juvenile"," Adult"),sep="")
    }
    if(length(SNnames)>1){
      xx<-rbind(SN_averages[[SNnames[1]]],SN_averages[[SNnames[2]]])
      xx<-xx[xx$SN>0,]
      thisSum<-summary(xx$SN)
      # med<-summary(xx$SN)[["Median"]]
    }else{
      xx<-SN_averages[[SNnames]]
      xx<-xx[xx$SN>0,]
      thisSum<-summary(xx$SN)
      # med<-summary(xx$SN)[["Median"]]
    }
    this_out<-list(med=thisSum[["Median"]],min=thisSum[["Min."]],max=thisSum[["Max."]])
    return(this_out)
  }
  get_KLP<-function(code){
    xx<-grep(paste("KLP_",code,sep=""),biolprm)
    thisKLP<-get_first_number(biolprm[xx])
    return(thisKLP)
  }
  get_KUP<-function(code){
    xx<-grep(paste("KUP_",code,sep=""),biolprm)
    thisKUP<-get_first_number(biolprm[xx])
    return(thisKUP)
  }

  if(sum(names(fun_group) == "InvertType") > 0)
    names(fun_group)[names(fun_group) == "InvertType"] <- "GroupType"
  
  # Subset invertebrates and vertebrates
  rs_names <- str_trim(fun_group[fun_group$GroupType %in% c("FISH", "MAMMAL", "SHARK", "BIRD"), "Name"]) # trim white space
  rs_codes <- str_trim(fun_group[fun_group$GroupType %in% c("FISH", "MAMMAL", "SHARK", "BIRD"), "Code"])
  invert_names <- fun_group[!(fun_group$GroupType %in% c("FISH", "MAMMAL", "SHARK", "BIRD")),]
  invert_codes<- str_trim(fun_group[!(fun_group$GroupType %in% c("FISH", "MAMMAL", "SHARK", "BIRD")), "Code"])

  cat("### ------------ Creating dynamic labels for vat                         ------------ ###\n")
  max_tracer <- nc_out$nvars
  max_layers <- length(nc_out$dim$z$vals)
  max_time <- length(nc_out$dim$t$vals)
  var_names <- names(nc_out$var)
  
  cat("### ------------ Creating map from BGM file                              ------------ ###\n")
  
  # Find number of boxes
  numboxes <- length(grep("# Box number", bgm))
  

  cat("### ------------ Setting up disaggregated spatial plots                  ------------ ###\n")
  nums <- grep("Nums", var_names, value = TRUE)
  N <- grep("_N", var_names, value = TRUE)
  N <- N[-grep("_Nums", N, value = FALSE)]
  tot_num <- c(nums)
  
  # extract tracers from the ncdf4 object
  vars <- list()
  for (i in 1:length(tot_num)){
    vars[[i]] <- ncvar_get(nc = nc_out, varid = tot_num[i])
  }
  names(vars) <- tot_num
  
  ###############
  ## RESN
  resns <- grep("ResN", var_names, value = TRUE)
  # extract tracers from the ncdf4 object

  RESN <- list()
  for (i in 1:length(resns)){
    RESN[[i]] <- ncvar_get(nc = nc_out, varid = resns[i])
  }
  names(RESN) <- resns
  ##################
  ## STRUCTN
  structns <- grep("StructN", var_names, value = TRUE)
  
  # extract tracers from the ncdf4 object
  STRUCTN <- list()
  for (i in 1:length(resns)){
    STRUCTN[[i]] <- ncvar_get(nc = nc_out, varid = structns[i])
  }
  names(STRUCTN) <- structns
  ##################
  
  
  # Create Erla's plots
  nominal_dz <- ncvar_get(nc = nc_out, varid = "nominal_dz")
  depth_layers <- nominal_dz[,which.max(colSums(nominal_dz))]
  depth_layers <- depth_layers[-c(1, length(depth_layers))]
  depth_layers <- cumsum(rev(depth_layers))
  
  depth_labels <- rep(NA, (length(depth_layers) + 1))
  for(i in 1:(length(depth_layers) + 1)){
    if(i == 1){
      depth_labels[i] <- paste("0 - ", depth_layers[i], sep = "")
    } else if(i == (length(depth_layers) + 1))
    {
      depth_labels[i] <- paste(depth_layers[i - 1], " + ", sep ="")
    } else depth_labels[i] <- paste(depth_layers[i - 1], " - ",  depth_layers[i], sep ="")
  }
  depth_labels <- c(depth_labels, "Sediment")
  
  vert_names <- fun_group[fun_group$GroupType %in% c("FISH", "MAMMAL", "SHARK", "BIRD"), "Code"]
  mat_age <- grep("_age_mat", biolprm, value = T)
  species_ids <- str_split_fixed(mat_age, "_age_mat", n = 2)[,1]
  #   juvenile_age <- as.numeric(gsub("[^\\d]+", "", mat_age, perl=TRUE))[which(species_ids %in% vert_names)]
  #test how many numbers there are. Sometimes there is an example value as well as the actual value, 
  #in which case we only want the first one
  source(paste(DIR$'General functions',"get_first_number.R",sep=""))
  # get_first_number<-function(x){
  #   yy<-gsub("([^\\d])","#",x,perl=TRUE)
  #   yyy<-unlist(str_split(yy,"#"))
  #   xPos<-grep("[^\\d]",yyy)[1]
  #   thisNum<-as.numeric(yyy[xPos])
  # }
  temp<-lapply(mat_age,FUN=get_first_number)
  juvenile_age<-as.numeric(unlist(temp))+1
  species_ids <- species_ids[which(species_ids %in% vert_names)]
  
  erla_plots <- list()
  SN_averages<-list()
  for(i in 1:length(species_ids)){
    spp <- fun_group[fun_group$Code == species_ids[i],c("Name", "NumCohorts")]
    spp <- str_trim(spp)
    if(juvenile_age[i] != 1){
      juv <- paste(spp[[1]], 1:(juvenile_age[i] - 1), "_Nums", sep = "")
      ad <- paste(spp[[1]], juvenile_age[i]:spp[[2]], "_Nums", sep = "")
      # Create the juveniles data
      juv_tmp <- NULL
      for(j in juv){
        x <- adply(vars[[j]], c(1, 3))
        juv_tmp <- rbind(juv_tmp, x)
      }
      juv_tmp <- juv_tmp %>%
        group_by(X1, X2) %>%
        summarize_each(funs(sum))
      colnames(juv_tmp) <- c("Layer", "Time", paste("Box", 0:(ncol(juv_tmp)-3), sep =" "))
      juv_tmp$Layer <- factor(juv_tmp$Layer,levels(juv_tmp$Layer)[c(((length(unique(juv_tmp$Layer)))-1):1, length(unique(juv_tmp$Layer)))])
      levels(juv_tmp$Layer) <- depth_labels  
      juv_tmp<-data.frame(juv_tmp) #turn to df so gather works
      juv_tmp <- gather(juv_tmp, Box, value = number, 3:ncol(juv_tmp))
      
      erla_plots[[paste(spp[[1]], "Juvenile")]] <- juv_tmp
      
      # Create the adults data
      ad_tmp <- NULL
      for(j in ad){
        x <- adply(vars[[j]], c(1, 3))
        ad_tmp <- rbind(ad_tmp, x)
      }
      
      ad_tmp <- ad_tmp %>%
        group_by(X1, X2) %>%
        summarize_each(funs(sum))
      colnames(ad_tmp) <- c("Layer", "Time", paste("Box", 0:(ncol(ad_tmp)-3), sep =" "))
      ad_tmp$Layer <- factor(ad_tmp$Layer,levels(ad_tmp$Layer)[c(((length(unique(ad_tmp$Layer)))-1):1, length(unique(ad_tmp$Layer)))])
      levels(ad_tmp$Layer) <- depth_labels
      ad_tmp<-data.frame(ad_tmp)
      ad_tmp <- gather(ad_tmp, Box, value = number, 3:ncol(ad_tmp))
      
      erla_plots[[paste(spp[[1]], "Adult")]] <- ad_tmp
      ###########################
      ## similar but for SN averages
      juv <- paste(spp[[1]], 1:(juvenile_age[i] - 1), "_StructN", sep = "")
      ad <- paste(spp[[1]], juvenile_age[i]:spp[[2]], "_StructN", sep = "")
      # Create the juveniles data
      juv_tmp <- NULL
      for(j in juv){
        x <- adply(STRUCTN[[j]], c(1, 3))
        juv_tmp <- rbind(juv_tmp, x)
      }
      juv_tmp <- juv_tmp %>%
        group_by(X1, X2) %>%
        summarize_each(funs(mean))
      colnames(juv_tmp) <- c("Layer", "Time", paste("Box", 0:(ncol(juv_tmp)-3), sep =" "))
      juv_tmp$Layer <- factor(juv_tmp$Layer,levels(juv_tmp$Layer)[c(((length(unique(juv_tmp$Layer)))-1):1, length(unique(juv_tmp$Layer)))])
      levels(juv_tmp$Layer) <- depth_labels  
      juv_tmp<-data.frame(juv_tmp) #turn to df so gather works
      juv_tmp <- gather(juv_tmp, Box, value = SN, 3:ncol(juv_tmp))
      
      SN_averages[[paste(spp[[1]], "Juvenile")]] <- juv_tmp
      
      # Create the adults data
      ad_tmp <- NULL
      for(j in ad){
        x <- adply(STRUCTN[[j]], c(1, 3))
        ad_tmp <- rbind(ad_tmp, x)
      }
      
      ad_tmp <- ad_tmp %>%
        group_by(X1, X2) %>%
        summarize_each(funs(mean))
      colnames(ad_tmp) <- c("Layer", "Time", paste("Box", 0:(ncol(ad_tmp)-3), sep =" "))
      ad_tmp$Layer <- factor(ad_tmp$Layer,levels(ad_tmp$Layer)[c(((length(unique(ad_tmp$Layer)))-1):1, length(unique(ad_tmp$Layer)))])
      levels(ad_tmp$Layer) <- depth_labels
      ad_tmp<-data.frame(ad_tmp)
      ad_tmp <- gather(ad_tmp, Box, value = SN, 3:ncol(ad_tmp))
      
      SN_averages[[paste(spp[[1]], "Adult")]] <- ad_tmp
      ####
      ############################################
    } 
    else {
      ad <- paste(spp[[1]], juvenile_age[i]:spp[[2]], "_Nums", sep = "")
      ad_tmp <- NULL
      for(j in ad){
        x <- adply(vars[[j]], c(1, 3))
        ad_tmp <- rbind(ad_tmp, x)
      }
      ad_tmp <- ad_tmp %>%
        group_by(X1, X2) %>%
        summarize_each(funs(sum))
      colnames(ad_tmp) <- c("Layer", "Time", paste("Box", 0:(ncol(ad_tmp)-3), sep =" "))
      ad_tmp$Layer <- factor(ad_tmp$Layer,levels(ad_tmp$Layer)[c(((length(unique(ad_tmp$Layer)))-1):1, length(unique(ad_tmp$Layer)))])
      levels(ad_tmp$Layer) <- depth_labels
      ad_tmp<-data.frame(ad_tmp)
      ad_tmp <- gather(ad_tmp, Box, value = number, 3:ncol(ad_tmp))
      erla_plots[[paste(spp[[1]], "Adult")]] <- ad_tmp
      ##################################
      ## do SN averages
      ad <- paste(spp[[1]], juvenile_age[i]:spp[[2]], "_StructN", sep = "")
      ad_tmp <- NULL
      for(j in ad){
        x <- adply(STRUCTN[[j]], c(1, 3))
        ad_tmp <- rbind(ad_tmp, x)
      }
      ad_tmp <- ad_tmp %>%
        group_by(X1, X2) %>%
        summarize_each(funs(mean))
      colnames(ad_tmp) <- c("Layer", "Time", paste("Box", 0:(ncol(ad_tmp)-3), sep =" "))
      ad_tmp$Layer <- factor(ad_tmp$Layer,levels(ad_tmp$Layer)[c(((length(unique(ad_tmp$Layer)))-1):1, length(unique(ad_tmp$Layer)))])
      levels(ad_tmp$Layer) <- depth_labels
      ad_tmp<-data.frame(ad_tmp)
      ad_tmp <- gather(ad_tmp, Box, value = number, 3:ncol(ad_tmp))
      SN_averages[[paste(spp[[1]], "Adult")]] <- ad_tmp
      ###
      #################################
    }
  }
  
  # --- End Erla Plots -- #
  
  age_species<-fun_group$Code[fun_group$NumCohorts>1]
  sngl_species<-fun_group$Code[fun_group$NumCohorts==1]
  all_species<-sort(c(age_species,sngl_species))
  
  pPREY <- data.frame(matrix(NA,ncol=(2*length(age_species)+length(sngl_species)),nrow=dim(fun_group)[1]))
  rownames(pPREY)<-fun_group$Code
  colnames(pPREY)<-sort(c(paste(sort(rep(age_species,2)),seq(1,2),sep=""),sngl_species))
  
  for(i in 1:length(all_species)){
    thisCode<-all_species[i]
    spp <- fun_group[fun_group$Code == thisCode,c("Name", "NumCohorts")]
    spp <- str_trim(spp)
    if(thisCode %in% age_species){
      for(age in 1:2){
        thisPreyVarNames<-paste("pPREY",seq(1,2),thisCode,age,sep="")
        xx<-str_trim(as.character(biolprm[grep(thisPreyVarNames[1],biolprm)+1]))
        xx<-as.double(unlist(str_split(xx," ")))
        yy<-str_trim(biolprm[grep(thisPreyVarNames[2],biolprm)+1])
        yy<-as.double(unlist(str_split(yy," ")))
        thisPreyAvail<-xx+yy
        thisColName<-paste(thisCode,age,sep="")
        pPREY[,match(thisColName,colnames(pPREY))]<-thisPreyAvail[1:dim(fun_group)[1]]
      }
    } else{
      thisPreyVarNames<-paste("pPREY",thisCode,sep="")
      testIndex<-grep(thisPreyVarNames,biolprm)
      if(length(testIndex)){
        xx<-str_trim(as.character(biolprm[testIndex+1]))
        xx<-as.double(unlist(str_split(xx," ")))
        thisPreyAvail<-xx
        thisColName<-paste(thisCode,sep="")
        pPREY[,match(thisColName,colnames(pPREY))]<-thisPreyAvail[1:dim(fun_group)[1]]
      }else{
        thisColName<-paste(thisCode,sep="")
        pPREY[,match(thisColName,colnames(pPREY))]<-rep(NA,dim(pPREY)[1])
      }
    }
  }
  
  # extract physical tracers from the ncdf4 object
  phy_names <- names(nc_out$var)[!(names(nc_out$var) %in% tot_num)]
  phy_names <- phy_names[-grep("_ResN", phy_names)] 
  phy_names <- phy_names[-grep("_StructN", phy_names)]
  vert_N <- paste(str_trim(rs_names), "_N", sep = "")
  phy_names <-  phy_names[!(phy_names %in% vert_N)]
  phy_names <- phy_names[-which(phy_names == "nominal_dz")]
  invert_nums <- grep("_N", phy_names, value = F)
  invert_mnames <- phy_names[invert_nums]
  trace_names <- phy_names[-(invert_nums)]
  
 
  invert_vars <- list()
  for (i in 1:length(invert_mnames)){
    tmp <- ncvar_get(nc = nc_out, varid = invert_mnames[i])
    if(length(dim(tmp)) == 2){
      if(dim(tmp)[1] == numboxes){
        tmp_invert <- tmp
        tmp_invert <- as.data.frame(tmp_invert)
        tmp_invert$Box <- paste("Box", 0:(numboxes - 1))
        tmp_invert <- gather(tmp_invert, Time, value = number, 1:(ncol(tmp_invert)-1))
        #keep just the number from Time
        tmp_invert$Time<-unlist(lapply(tmp_invert$Time,get_first_number))
        levels(tmp_invert$Time) <- 0:length(unique(tmp_invert$Time))
        tmp_invert$Time <- as.numeric(as.character(tmp_invert$Time))
        # refErla<-seq(1,(length(erla_plots)))[names(erla_plots)==invert_mnames[i]]
        # erla_plots[[refErla]] <- tmp_invert
        
        erla_plots[[invert_mnames[i]]] <- tmp_invert
        invert_vars[[i]] <- tmp
      }
    } else{
      tmp_invert <- adply(tmp, c(1,3))
      tmp_invert <- tmp_invert %>%
        group_by(X1, X2) %>%
        summarize_each(funs(sum))
 
      colnames(tmp_invert) <- c("Layer", "Time", paste("Box", 0:(ncol(tmp_invert)-3), sep =" "))
      tmp_invert$Layer <- factor(tmp_invert$Layer,levels(tmp_invert$Layer)[c(((length(unique(tmp_invert$Layer)))-1):1, length(unique(tmp_invert$Layer)))])
      levels(tmp_invert$Layer) <- depth_labels
      tmp_invert<-data.frame(tmp_invert)
      tmp_invert <- gather(tmp_invert, Box, value = number, 3:ncol(tmp_invert))
      invert_vars[[i]] <- tmp
      # refErla<-seq(1,(length(erla_plots)))[names(erla_plots)==invert_mnames[i]]
      # if(length(refErla)>0){
      #   erla_plots[[refErla]] <- tmp_invert
      # }
      erla_plots[[invert_mnames[i]]] <- tmp_invert
    }
  }
  names(invert_vars) <- invert_mnames
  
  trace_vars <- list()
  for (i in 1:length(trace_names)){
    tmp <- ncvar_get(nc = nc_out, varid = trace_names[i])
    if(length(dim(tmp)) == 2){
      if(dim(tmp)[1] == numboxes){
        tmp_trace <- tmp
        tmp_trace <- as.data.frame(tmp_trace)
        tmp_trace$Box <- paste("Box", 0:(numboxes - 1))
        tmp_trace<-data.frame(tmp_trace)
        tmp_trace <- gather(tmp_trace, Time, value = number, 1:(ncol(tmp_trace)-1))
        levels(tmp_trace$Time) <- 0:length(unique(tmp_trace$Time))
        tmp_trace$Time <- as.numeric(as.character(tmp_trace$Time))
        erla_plots[[trace_names[i]]] <- tmp_trace
        trace_vars[[i]] <- tmp
      }
    } else{
      tmp_array <- adply(tmp, c(1,3))
      tmp_trace <- tmp_array %>%
        group_by(X1, X2) %>%
        summarize_each(funs(sum))
      colnames(tmp_trace) <- c("Layer", "Time", paste("Box", 0:(ncol(tmp_trace)-3), sep =" "))
      tmp_trace$Layer <- factor(tmp_trace$Layer,levels(tmp_trace$Layer)[c(((length(unique(tmp_trace$Layer)))-1):1, length(unique(tmp_trace$Layer)))])
      levels(tmp_trace$Layer) <- depth_labels
      tmp_trace<-data.frame(tmp_trace)
      tmp_trace <- gather(tmp_trace, Box, number, 3:ncol(tmp_trace))
      
      trace_vars[[i]] <- tmp
      erla_plots[[trace_names[i]]] <- tmp_trace
    }
  }
  names(trace_vars) <- trace_names
  
  
##########################
  ## predator-prey encounter rates
  orderBox<-data.frame(cbind("seq"<-seq(1,numboxes),"box"<-paste("Box.",seq(0,(numboxes-1)),sep="")))
  
    PREY_SELECT<-NULL
    PRED_SELECT<-NULL
    
    PREY_SELECT_SN<-NULL
    PRED_SELECT_SN<-NULL
    
    PREY_SELECT_PREYAVAIL<-NULL
    PRED_SELECT_PREYAVAIL<-NULL

    fgcodes<-fun_group[,"Code"]
    fgnames <- fun_group[,"Name"]

    for(g in 1:(length(fgcodes))){
      thisCode<-fgcodes[g]
      thisName<-str_trim(fgnames[g],side="both")
      cat("\n",thisName,"-->")
      if(thisCode %in% invert_names$Code){
        thisNumCohorts<-invert_names$NumCohorts[invert_names$Code==thisCode]
        if(thisNumCohorts>1){
          thisVars<-paste(thisName,"_N",seq(1,2),sep="")
          #add up the data for both juvenile and adult 
          for(j in 1:(length(thisVars))){
            refErla<-seq(1,(length(erla_plots)))[names(erla_plots)==thisVars[j]]
            thisData<-erla_plots[[refErla]]
            
            #check if we need to scale up the dimensions to include layers
            if(dim(thisData)[1] == max_time * numboxes){
              combLayerData<-data.frame(matrix(NA,ncol=4,nrow=0))
              colnames(combLayerData)<-c("Layer","Box","Time","number")
              for(l in 1:(length(depth_labels))){
                thisLayer<-depth_labels[l]
                thisLayerData<-cbind("Layer"=rep(thisLayer,dim(thisData)[1]),thisData)
                if(l!=max_layers){
                  thisLayerData$number<-0*thisLayerData$number
                }
                combLayerData<-rbind(combLayerData,thisLayerData)
              }
              thisData<-combLayerData
              rm(combLayerData)
              
            }
            
            if(j==1){
              tmp <- thisData
 
            } else{
              tmp$number<-tmp$number+thisData$number
            }
          }
        } else{
          thisVar<-paste(thisName,"_N",sep="")
          refErla<-seq(1,(length(erla_plots)))[names(erla_plots)==thisVar]
          
          thisData<-erla_plots[[refErla]]
          #check here if dim(thisData)[1] == maxlayers * numBoxes * numtimesteps
          if(dim(thisData)[1] == max_time * numboxes){
            combLayerData<-data.frame(matrix(NA,ncol=4,nrow=0))
            colnames(combLayerData)<-c("Layer","Box","Time","number")
            for(l in 1:(length(depth_labels))){
              thisLayer<-depth_labels[l]
              thisLayerData<-cbind("Layer"=rep(thisLayer,dim(thisData)[1]),thisData)
              if(l!=max_layers){
                thisLayerData$number<-0*thisLayerData$number
              }
              combLayerData<-rbind(combLayerData,thisLayerData)
            }
            thisData<-combLayerData
            rm(combLayerData)

          }
          tmp <- thisData
        }
        
      } else{
        thisVars<-paste(thisName,c(" Juvenile"," Adult"),sep="")
        #add up the data for both juvenile and adult 
        for(j in 1:(length(thisVars))){
          refErla<-seq(1,(length(erla_plots)))[names(erla_plots)==thisVars[j]]
          thisData<-erla_plots[[refErla]]
          if(j==1){
            tmp <- thisData
          } else{
            tmp$number<-tmp$number+thisData$number
          }
        }
      }
      
      tmp$Time <- as.numeric(as.character(tmp$Time)) * toutinc / 365 + startyear
      
      ##DO THE PREDATORS OF THIS GROUP FIRST
      cat("Getting predators...")
      #this is the potential predators for the selected prey
      pPRED<-pPREY[rownames(pPREY)==thisCode,]
      index<-pPRED>0 & !is.na(pPRED)
      if(length(pPRED[index])>0){
        thisPredators<-colnames(pPRED)[index]
        thisPredatorNames<-thisPredators
        
        thisPreyAvail<-data.frame(cbind(colnames(pPRED)[index],pPRED[index]))
        colnames(thisPreyAvail)<-c("Predator","Available")
        PREY_SELECT_PREYAVAIL[[thisCode]]<-thisPreyAvail
        
        #get SN for predators and prey and record average ratio
        #check if the prey is a vertebrate first. If it isn't it doesn't have SN
        if(thisCode %in% vert_names){
          thisSNsum<-getAverageSN(code=thisCode,age="both")
          thisSN<-thisSNsum$'med'
          ppSN<-data.frame(matrix(NA,ncol=5,nrow=length(thisPredators)))
          colnames(ppSN)<-c("Predator","SN","SNratio","LowerGape","UpperGape")
          ppSN$Predator<-thisPredators
          for(i in 1:(length(thisPredators))){
            xx<-get_first_number(thisPredators[i])
            predAge<-ifelse(xx==2,"Adult","Juvenile")
            predCode<-unlist(str_split(thisPredators[i],as.character(xx)))[1]
            #check if predator is vertebrate.
            #if its not, need to get sn from biolprm
            if(predCode %in% vert_names){
              predSNsum<-getAverageSN(predCode,predAge)
              predSN<-predSNsum$'med'
              predSNmax<-predSNsum$'max'
              predSNmin<-predSNsum$'min'
            }else{
              xx<-grep(paste("^",predCode,"_sn",sep=""),biolprm)
              predSN<-get_first_number(biolprm[xx])
              predSNmax<-predSN; predSNmin<-predSN;
            }
            ppSN$SN[i]<-predSN
            ppSN$SNratio[i]<-thisSN/predSN
            predKLP<-get_KLP(predCode)
            predKUP<-get_KUP(predCode)
            ppSN$LowerGape[i]<-predKLP*predSN
            ppSN$UpperGape[i]<-predKUP*predSN
            
            ppSN$LowerGapemin[i]<-predKLP*predSNmin
            ppSN$UpperGapemax[i]<-predKUP*predSNmax
          }
          ppSN$preySN<-thisSN
          ppSN$preySN_min<-thisSNsum$'min'
          ppSN$preySN_max<-thisSNsum$'max'
          PREY_SELECT_SN[[thisCode]]<-ppSN
        }
        
        #get data for all other species and testing for overlap
        test<-tmp
        
        test$FullInt<-rep(1,nrow(test))
        fullByBox<-tapply(test$FullInt,test$Box,sum,na.rm=TRUE)
        
        byBox<-data.frame(matrix(NA,nrow=numboxes,ncol=length(thisPredators)))
        colnames(byBox)<-thisPredators
        
        for(i in 1:(length(thisPredators))){
          thisPredCode<-gsub("[\\d]","",thisPredators[i],perl=TRUE)
          ##NEED TO HANDLE INVERTS IN HERE
          #######################
          if(!(thisPredCode %in% age_species)){
            thisPredName<-str_trim(fgnames[match(thisPredCode,fgcodes)])
            varPredName<-paste(thisPredName,"_N",sep="")
          } else{
            thisPredAge<-get_first_number(thisPredators[i])
            thisPredName<-str_trim(fgnames[match(thisPredCode,fgcodes)])
            if(thisPredCode %in% invert_codes){
              varPredName<-paste(thisPredName,"_N",thisPredAge,sep="")
            } else {
              varPredName<-ifelse(thisPredAge==1,paste(thisPredName," Juvenile",sep=""),paste(thisPredName," Adult",sep=""))
            }
          }
          cat(varPredName,"-->")
          thisI<-seq(1,(length(erla_plots)))[names(erla_plots)==varPredName]
          thisData<-erla_plots[[thisI]]
          
          if(dim(thisData)[1]==dim(test)[1]){
            
            test$numberPred<-thisData$number
            test$interact<-mapply(FUN=function(x,y){ifelse((x>0) & (y>0),1,0)},tmp$number,thisData$number)
            
            thisByBox<-tapply(test$interact,test$Box,sum,na.rm=TRUE)
            byBox[,i]<-thisByBox/fullByBox
          } else if(dim(thisData)[1]*max_layers==dim(test)[1]){
            thisData$Time <- as.numeric(as.character(thisData$Time)) * toutinc / 365 + startyear
            test$number<-rep(NA,dim(test)[1])
            for(t in unique(test$Time)){
              for(b in unique(test$Box)){
                testIndex<-test$Box==b & test$Time==t
                thisIndex<-thisData$Box==b & thisData$Time==t
                if(length(thisData$number[thisIndex])>0){
                  test$number[testIndex]<-thisData$number[thisIndex]
                }
              }
            }
            test$interact<-mapply(FUN=function(x,y){ifelse((x>0) & (y>0),1,0)},tmp$number,thisData$number)
            
            thisByBox<-tapply(test$interact,test$Box,sum,na.rm=TRUE)
            byBox[,i]<-thisByBox/fullByBox
          } else{
            byBox[,i]<-0*byBox[,i]
          }
          
          
        }
        PREY_SELECT[[thisCode]]<-byBox
      }
      
      ##NOW DO THE PREYS OF THIS GROUP FIRST
      cat("Getting preys ...")
      #this is the potential predators for the selected prey
      pPRED<-pPREY[,(colnames(pPREY) %in% paste(thisCode,c("","1","2"),sep=""))] ##NEED TO ADD INVERTS INTO obj$pPREY AS PREDATORS
      index<-!is.na(pPRED)
      test<-pPRED[index]
      if(length(test)>0){
        if(length(dim(pPRED))>0){
          pPRED$total<-apply(pPRED,1,sum)
        } else{
          pPRED<-data.frame(cbind(thisCode=pPRED,"total"=pPRED))
          colnames(pPRED)<-c(thisCode,"total")
          rownames(pPRED)<-rownames(pPREY)
        }
        index<-pPRED$total>0 & !is.na(pPRED$total)
        if(length(pPRED[index,1])>0){
          thisPreys<-rownames(pPRED)[index]
          thisPreyNames<-thisPreys
          
          #get SN for predator and preys and record average ratio
          #check if the predator is a vertebrate first. If it isn't it doesn't have SN
          if(thisCode %in% vert_names){
            thisSNsum<-getAverageSN(code=thisCode,age="both")
            thisSN<-thisSNsum$'med'
            thisSNmin<-thisSNsum$'min'
            thisSNmax<-thisSNsum$'max'
          }else{
            xx<-grep(paste("^",thisCode,"_sn",sep=""),biolprm)
            thisSN<-get_first_number(biolprm[xx])
            thisSNmin<-thisSN
            thisSNmax<-thisSN
          }
            ppSN<-data.frame(matrix(NA,ncol=5,nrow=length(thisPreys)))
            colnames(ppSN)<-c("Prey","SN","SNratio","LowerGape","UpperGape")
            ppSN$Prey<-thisPreys
            predKLP<-get_KLP(thisCode)
            predKUP<-get_KUP(thisCode)
            ppSN$LowerGape<-predKLP*thisSN
            ppSN$UpperGape<-predKUP*thisSN
            
            ppSN$UpperGapemax<-predKUP*thisSNmax
            ppSN$LowerGapemin<-predKLP*thisSNmin
            
            for(i in 1:(length(thisPreys))){
              xx<-get_first_number(thisPreys[i])
              preyAge<-ifelse(is.na(xx),"both",ifelse(xx==2,"Adult","Juvenile"))
              if(is.na(xx)){
                preyCode<-thisPreys[i]
              }else{
                preyCode<-unlist(str_split(thisPreys[i],as.character(xx)))[1]
              }
              #check if prey is an invert - if it is not it doesn't matter what gapesize is
              if(preyCode %in% vert_names){
                preySNsum<-getAverageSN(preyCode,preyAge)
                preySN<-preySNsum$'med'
                preySNmin<-preySNsum$'min'
                preySNmax<-preySNsum$'max'
              } else{
                preySN<-NA
                preySNmin<-NA
                preySNmax<-NA
              }
              ppSN$SN[i]<-preySN
              ppSN$SNmin[i]<-preySNmin
              ppSN$SNmax[i]<-preySNmax
              ppSN$SNratio[i]<-preySN/thisSN

            }
            ppSN$predSN<-thisSN
            ppSN$predSNmin<-thisSNmin
            ppSN$predSNmax<-thisSNmax
            PRED_SELECT_SN[[thisCode]]<-ppSN
        } else{
            thisPreys<-c()
        }
          ##########################################
          
          ###edit this
          thisPreyAvail<-data.frame(cbind(thisPreys,pPRED[index,1]))
          colnames(thisPreyAvail)<-c("Prey","Available")
          PRED_SELECT_PREYAVAIL[[thisCode]]<-thisPreyAvail
          
          #get data for all available prey species and test for overlap
          test<-tmp
          
          test$FullInt<-rep(1,nrow(test))
          fullByBox<-tapply(test$FullInt,test$Box,sum,na.rm=TRUE)
          
          byBox<-data.frame(matrix(NA,nrow=numboxes,ncol=length(thisPreys)))
          colnames(byBox)<-thisPreys
          
          if(length(thisPreys)>0){
            for(i in 1:(length(thisPreys))){
            thisPreyCode<-gsub("[\\d]","",thisPreys[i],perl=TRUE)
            thisPreyName<-str_trim(fgnames[match(thisPreyCode,fgcodes)])
            cat(thisPreyName,"-->")
            #I have prey grouped together for juveniles and adults, so will need to check to encounters with both where the prey is a vertebrate
            #DO INVERT PREYS FIRST
            if(thisPreyCode %in% invert_names$Code){
              #WITH AGE-STRUCTURE
              thisNumCohorts<-invert_names$NumCohorts[invert_names$Code==thisPreyCode]
              if(thisNumCohorts>1){
                test$numberPrey1<-0*(test$number)
                test$numberPrey2<-0*(test$number)
                for(age in 1:2){
                  varPreyName<-ifelse(age==1,paste(thisPreyName,"_N",age,sep=""),paste(thisPreyName,"_N",age,sep=""))
                  thisI<-seq(1,(length(erla_plots)))[names(erla_plots)==varPreyName]
                  thisData<-erla_plots[[thisI]]
                  
                  #  if(dim(thisData)[1]==dim(test)[1]){
                  
                  test[,paste("numberPrey",age,sep="")]<-thisData$number
                  #  }
                  
                }
                #this bit is slow. replace..?
                test$interact<-mapply(FUN=function(x,y,z){ifelse((x>0) & ((y>0) | (z>0)),1,0)},test$number,test$numberPrey1,test$numberPrey2)
                
                thisByBox<-tapply(test$interact,test$Box,sum,na.rm=TRUE)

                thisByBoxFullBox<-thisByBox/fullByBox
                #reorder
                names(orderBox)<-c("seq","box")
                orderThisBox<-thisByBoxFullBox[match(orderBox$box,names(thisByBoxFullBox))]
                
                byBox[,i]<-orderThisBox
                
              } else{
                #INVERT PREYS WITHOUT AGE-STRUCTURE
                varPreyName<-paste(thisPreyName,"_N",sep="")
                thisI<-seq(1,(length(erla_plots)))[names(erla_plots)==varPreyName]
                thisData<-erla_plots[[thisI]]
                
                if(dim(thisData)[1]==dim(test)[1]){
                  
                  test[,"numberPrey"]<-thisData$number
                  
                  test$interact<-mapply(FUN=function(x,y){ifelse((x>0) & (y>0),1,0)},test$number,test$numberPrey)
                  
                  thisByBox<-tapply(test$interact,test$Box,sum,na.rm=TRUE)

                  thisByBoxFullBox<-thisByBox/fullByBox
                  #reorder
                  names(orderBox)<-c("seq","box")
                  orderThisBox<-thisByBoxFullBox[match(orderBox$box,names(thisByBoxFullBox))]
                  
                  byBox[,i]<-orderThisBox
                  
                } else if(dim(thisData)[1]*max_layers==dim(test)[1]){
                  thisData$Time <- as.numeric(as.character(thisData$Time)) * toutinc / 365 + startyear
                  test$numberPrey<-rep(NA,dim(test)[1])
                  for(t in unique(test$Time)){
                    for(b in unique(test$Box)){
                      altb<-gsub("\\."," ",b,perl=TRUE)
                      testIndex<-(test$Box==b | test$Box==altb) & test$Time==t
                      thisIndex<-(thisData$Box==b | thisData$Box==altb) & thisData$Time==t
                      if(length(thisData$number[thisIndex])>0){
                        test$numberPrey[testIndex]<-thisData$number[thisIndex]
                      } else{
                        test$numberPrey[testIndex]<-rep(NA,length(test$number[testIndex]))
                      }
                    }
                  }
                  test<-data.frame(test)
                  test$interact<-mapply(FUN=function(x,y){ifelse((x>0) & (y>0),1,0)},test$numberPrey,test$number)
                  
                  thisByBox<-tapply(test$interact,test$Box,sum,na.rm=TRUE)
 
                  thisByBoxFullBox<-thisByBox/fullByBox
                  #reorder
                  names(orderBox)<-c("seq","box")
                  orderThisBox<-thisByBoxFullBox[match(orderBox$box,names(thisByBoxFullBox))]
                  
                  byBox[,i]<-orderThisBox
                  
                } 
                
              }
              
              
            } else{
              #NOW DO VERTEBRATE PREYS
              test$numberPrey1<-0*(test$number)
              test$numberPrey2<-0*(test$number)
              for(age in 1:2){
                varPreyName<-ifelse(age==1,paste(thisPreyName," Juvenile",sep=""),paste(thisPreyName," Adult",sep=""))
                thisI<-seq(1,(length(erla_plots)))[names(erla_plots)==varPreyName]
                thisData<-erla_plots[[thisI]]
                test[,paste("numberPrey",age,sep="")]<-thisData$number
              }
              
              test$interact<-mapply(FUN=function(x,y,z){ifelse((x>0) & ((y>0) | (z>0)),1,0)},test$number,test$numberPrey1,test$numberPrey2)
              
              thisByBox<-tapply(test$interact,test$Box,sum,na.rm=TRUE)
              thisByBoxFullBox<-thisByBox/fullByBox
              #reorder
              names(orderBox)<-c("seq","box")
              orderThisBox<-thisByBoxFullBox[match(orderBox$box,names(thisByBoxFullBox))]
              
              byBox[,i]<-orderThisBox
              
            }
            
            }
          }
          
          PRED_SELECT[[thisCode]]<-byBox 
        }
    }

  
  
  output <- list(disagg = vars,invert_vars = invert_vars, invert_mnames = invert_mnames, trace_vars = trace_vars, trace_names = trace_names, var_names = tot_num, max_layers = max_layers, max_time = max_time,  rs_names = rs_names, islands = islands, numboxes = numboxes, fun_group = fun_group, invert_names = invert_names,erla_plots = erla_plots, toutinc = toutinc, startyear = startyear, fgnames = fgnames,fgcodes=fgcodes,pPREY=pPREY,PREY_SELECT=PREY_SELECT,PRED_SELECT=PRED_SELECT,PREY_SELECT_SN=PREY_SELECT_SN,PRED_SELECT_SN=PRED_SELECT_SN,PREY_SELECT_PREYAVAIL=PREY_SELECT_PREYAVAIL,PRED_SELECT_PREYAVAIL=PRED_SELECT_PREYAVAIL)

  cat("### ------------ vat object created, you can now run the vat application ------------ ###\n") 
  return(output)
  class(output) <- "vadt"
}