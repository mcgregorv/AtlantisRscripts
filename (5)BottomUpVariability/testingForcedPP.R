#read in tracers and check out temperature to see if correct ROMS being used
this_run<-"base"
this_out<-"TestPP"
baseOut<-"BaseNEW"

basePath<-paste(DIR$'Base',"ATLANTISmodels\\",sep="")

ThisNC.nc<-nc_open(paste(basePath,"base\\output",this_out,"\\output.nc", sep=""))
BaseNC.nc<-nc_open(paste(basePath,"base\\",baseOut,"\\output.nc", sep=""))

thisVol<-ncvar_get(ThisNC.nc, "volume")
nlayers<-dim(thisVol)[1]
baseVol<-ncvar_get(BaseNC.nc, "volume")

thisDiatom<-ncvar_get(ThisNC.nc, "Diatom_N")
baseDiatom<-ncvar_get(BaseNC.nc, "Diatom_N")

dynBoxes<-2:25

thisDiatomBiomass<-apply(thisDiatom[-nlayers,dynBoxes,] * thisVol[-nlayers, dynBoxes,], 3, sum) * mg_2_tonne * X_CN
baseDiatomBiomass<-apply(baseDiatom[-nlayers, dynBoxes, ] * baseVol[-nlayers, dynBoxes, ], 3, sum) * mg_2_tonne * X_CN

thisMeanDiatom<-apply(thisDiatom[-nlayers, dynBoxes,], 3, nonZeroMean)
baseMeanDiatom<-apply(baseDiatom[-nlayers, dynBoxes,], 3, nonZeroMean)

allTracers<-names(BaseNC.nc$var)
