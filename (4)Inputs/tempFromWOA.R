##these data are from WOA https://www.nodc.noaa.gov/OC5/woa13/woa13data.html
#3 netcdf files for year ranges 1985-1994; 1995-2004; 2005-2012

dataPath<-paste(DIR$'Data',"\\ocean_physics\\temperature\\",sep="")

tempNC1.nc<-nc_open(paste(dataPath,"woa13_8594_t01_01v2.nc",sep=""))
tempNC2.nc<-nc_open(paste(dataPath,"woa13_95A4_t01_01v2.nc",sep=""))
tempNC3.nc<-nc_open(paste(dataPath,"woa13_A5B2_t01_01v2.nc",sep=""))

aveDec<-nc_open(paste(dataPath,"woa13_decav_t01_01v2.nc",sep=""))

vars<-names(tempNC1.nc$var)
# vars
# [1] "crs"                "lat_bnds"           "lon_bnds"           "depth_bnds"         "climatology_bounds" "t_an"              
# [7] "t_mn"               "t_dd"               "t_sd"               "t_se"               "t_oa"               "t_ma"              
# [13] "t_gp"
##definitions here https://data.nodc.noaa.gov/woa/WOA13/DOC/woa13documentation.pdf
##t_an objectively analysed sea water temperature
##t_mn is Average of all unflagged interpolated sea_water_temperature

##I'm going to use t_an. we also need lat, lon. There is no time here, it is just for January. Need depth too

temp1<-ncvar_get(tempNC1.nc,"t_an"); temp2<-ncvar_get(tempNC2.nc,"t_an"); temp3<-ncvar_get(tempNC3.nc,"t_an")
lat1<-ncvar_get(tempNC1.nc,"lat_bnds"); lat2<-ncvar_get(tempNC1.nc,"lat_bnds"); lat3<-ncvar_get(tempNC1.nc,"lat_bnds")
lon1<-ncvar_get(tempNC1.nc,"lon_bnds"); lon2<-ncvar_get(tempNC1.nc,"lon_bnds"); lon3<-ncvar_get(tempNC1.nc,"lon_bnds")
depth1<-ncvar_get(tempNC1.nc,"depth_bnds"); depth2<-ncvar_get(tempNC1.nc,"depth_bnds"); depth3<-ncvar_get(tempNC1.nc,"depth_bnds")

##check out the dimensions. lat, lon, and depth are 1D. temperature is 3D and its dimensions are the lengths of lat,lon,depth (respectively)
# dim(temp1)
# [1] 360 180  57
# > length(lat1)
# [1] 360
# > length(lon1)
# [1] 720
# > length(depth1)
# [1] 114


shapeFile<-paste(dataPath,"\\sss\\sss",sep="")

shape<-read.shapefile(shapeFile)






