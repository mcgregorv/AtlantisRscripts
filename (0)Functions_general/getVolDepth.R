getVolDepth<-function(bgmFile,depthLayers){
  
  #read in the bgm file
  bgm<-read_boxes(fnm=bgmFile)
  this_nbox<-bgm$'nbox'
  
  dl_sort<-sort(depthLayers)
  dl<-depthLayers
  num_wc<-length(dl)
  lt<-c(rep("wc",(length(dl))),rep("sd",1)) # one sediment, rest WC
  
  #for each box and for each depth bin, need to calculate the volume. We have the area of each box bgm$'b_area
  #also have the mean depth for each box (bgm$'b_botz')
  vol<-data.frame(matrix(NA,nrow=this_nbox,ncol=(length(dl)+3)))
  colnames(vol)<-c("varname","boxid",paste("dl",seq(0,(num_wc-1)),sep=""),"sl0")
  vol$varname<-"volume"
  vol$boxid<-seq(0,(this_nbox-1))
  dz<-vol
  dz$varname<-rep("dz",nrow(vol))
  for(b in 1:bgm$'nbox'){
    #   cat(paste(b,"\n"))
    this_area<-as.double(bgm$'b_area'[b])
    this_mean_depth<-abs(as.double(bgm$'b_botz'[b]))
    #which depth layers are in this box?
    if(this_mean_depth>max(dl)){
      this_dl<-c(this_mean_depth,sort(dl_sort,decreasing=TRUE)[-1],1)
      #now fix this_dl so that for the 2nd entry through to the 2nd to last entry are the difference between the current and the previous
      new_dl<-this_dl
      for(i in 1:(length(this_dl)-2)){
        new_dl[i]<-this_dl[i]-this_dl[i+1]
      }
      this_dl<-new_dl
      rm(new_dl)
    }else{
      temp<-dl_sort<this_mean_depth
      this_dl<-dl_sort[temp]
      #add the next depth layer with the remaining depth
      next_layer<-this_mean_depth-max(c(this_dl,0))
      this_dl<-c(next_layer,sort(this_dl,decreasing=TRUE))
    
      #then if need be, add 0's for any layers not used
      this_dl<-c(this_dl,rep(0,(num_wc-length(this_dl))),1)
      #now fix this_dl so that for the 2nd entry through to the 2nd to last entry are the difference between the current and the previous
      new_dl<-this_dl
      for(i in 2:(length(this_dl)-2)){
        new_dl[i]<-this_dl[i]-this_dl[i+1]
      }
      this_dl<-c(new_dl) 
      rm(new_dl)
    }
  
    dz[b,3:(ncol(dz))]<-this_dl
    for(d in 1:(length(this_dl))){
      this_depth<-this_dl[d]
      this_type<-lt[d]
      vol[b,(d+2)]<-this_depth*this_area
    }
  }
  
  nominal_dz<-dz
  nominal_dz$varname<-rep("nominal_dz",nrow(vol))
  
  layer_temp<-nominal_dz[,3:ncol(nominal_dz)]
  layer_index<-layer_temp>0
  rm(layer_temp)
  returnThese<-list("bgm"=bgm,"vol"=vol,"dz"=dz,"nominal_dz"=nominal_dz,"layer_index"=layer_index)
  return(returnThese)
}