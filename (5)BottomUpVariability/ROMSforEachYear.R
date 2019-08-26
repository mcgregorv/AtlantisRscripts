#have ROMS outputs for 9 years - 1996-2004. Create ROMS .nc files for each year so these can be used for seperate models to test effect of repeating any one year
#deifne base ROMS file which has one year in it - I think 2004
romsYears<-seq(1996,2004)

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")
version<-"A"

nlayer<-6; nboxes<-30

#read in the full ROMS files
ROMSpath<-paste(basePath,"inputs\\ROMS\\",sep=""); 
inputROMSPath<-ROMSpath
# inputROMSPath<-paste(ROMSpath,"BootstrapROMS\\",version,"\\",sep="")


temp1yearROMSFile<-paste(ROMSpath,"Chatham30_temp.nc",sep="")
hydro1yearROMSFile<-paste(ROMSpath,"Chatham30_hydro.nc",sep="")
salt1yearROMSFile<-paste(ROMSpath,"Chatham30_salt.nc",sep="")

tempFullROMSFile<-paste(ROMSpath,"Chatham30_tempAll.nc",sep="")
hydroFullROMSFile<-paste(ROMSpath,"Chatham30_hydroAll.nc",sep="")
saltFullROMSFile<-paste(ROMSpath,"Chatham30_saltAll.nc",sep="")

fullHydro<-nc_open(hydroFullROMSFile)
fullTemp<-nc_open(tempFullROMSFile)
fullSalt<-nc_open(saltFullROMSFile)

vnames_salt<-names(fullSalt$var); vnames_temp<-names(fullTemp$var); vnames_hydro<-names(fullHydro$var)
saltData<-ncvar_get(fullSalt,"salinity"); tempData<-ncvar_get(fullTemp,"temperature")
nts<-dim(saltData)[3] #they are 12 hour time steps
hydroExchangeData<-ncvar_get(fullHydro,"exchange"); hydroDestBData<-ncvar_get(fullHydro,"dest_b"); hydroDestKData<-ncvar_get(fullHydro,"dest_k")

nyears<-nts/(365*2) #check it is a round number!
tsPerYear<-2*365
#split the data up by year
saltByYear<-NULL; tempByYear<-NULL; exchangeByYear<-NULL; destBByYear<-NULL;  destKByYear<-NULL
for(y in 1:nyears){
  startY<-(y-1)*tsPerYear+1; endY<-startY+tsPerYear-1
  # cat(paste(y, " has start year",startY,", endY ", endY," which in years are ",startY/(365*2)," ",endY/(365*2),"\n",sep=""))
  saltByYear[[y]]<-saltData[,,startY:endY]; tempByYear[[y]]<-tempData[,,startY:endY]; 
  exchangeByYear[[y]]<-hydroExchangeData[,,,startY:endY]
  destBByYear[[y]]<-hydroDestBData[,,,startY:endY]; destKByYear[[y]]<-hydroDestKData[,,,startY:endY]
}

for(y in 1:nyears){
  thisYear<-romsYears[y]
  ##TEMPERATURE FILES
  newTempNCfile<-paste(inputROMSPath,"Chatham_temp_year",thisYear,".nc",sep="")
  file.copy(from=temp1yearROMSFile,to=newTempNCfile,overwrite=TRUE)
  newTempNC<-nc_open(newTempNCfile,write=TRUE)
  thisValsData<-tempByYear[[y]]
  ncvar_put(newTempNC,varid = "temperature",vals = thisValsData)
  nc_close(newTempNC)
  ##SALINITY FILES
  newSaltNCfile<-paste(inputROMSPath,"Chatham_salt_year",thisYear,".nc",sep="")
  file.copy(from=salt1yearROMSFile,to=newSaltNCfile,overwrite=TRUE)
  newSaltNC<-nc_open(newSaltNCfile,write=TRUE)
  thisValsData<-saltByYear[[y]]
  ncvar_put(newSaltNC,varid = "salinity",vals = thisValsData)
  nc_close(newSaltNC)
  ##HYDRO FILES vars "exchange" "dest_b"   "dest_k"  
  newHydroNCfile<-paste(inputROMSPath,"Chatham_hydro_year",thisYear,".nc",sep="")
  file.copy(from=hydro1yearROMSFile,to=newHydroNCfile,overwrite=TRUE)
  newHydroNC<-nc_open(newHydroNCfile,write=TRUE)
  thisValsData<-exchangeByYear[[y]]
  ncvar_put(newHydroNC,varid = "exchange",vals = thisValsData)
  ##
  thisValsData<-destBByYear[[y]]
  ncvar_put(newHydroNC,varid = "dest_b",vals = thisValsData)
  ##
  thisValsData<-destKByYear[[y]]
  ncvar_put(newHydroNC,varid = "dest_k",vals = thisValsData)
  #close the nc file
  nc_close(newHydroNC)
}

#####################################
## TESTING ###############################
# read in hydro all years file, and compare with individual hydro files
fullExchange<-ncvar_get(fullHydro,"exchange")
fullDestb<-ncvar_get(fullHydro,"dest_b"); fullDestk<-ncvar_get(fullHydro,"dest_k")
yearTimeSteps<-sort(rep(romsYears,730))
y<-1; thisYear<-romsYears[y]; yTimestepIndex<-yearTimeSteps == thisYear
thisExchange<-fullExchange[,,,yTimestepIndex]
thisDestb<-fullDestb[,,,yTimestepIndex]; thisDestk<-fullDestk[,,,yTimestepIndex]

yearExchange<-nc_open(paste(ROMSpath,"Chatham_hydro_year",thisYear,".nc",sep=""))
yearHydro<-ncvar_get(yearExchange,"exchange")
yearDestb<-ncvar_get(yearExchange,"dest_b"); yearDestk<-ncvar_get(yearExchange,"dest_k")

all.equal(yearHydro, thisExchange)
all.equal(yearDestb, thisDestb)
all.equal(yearDestk, thisDestk)
