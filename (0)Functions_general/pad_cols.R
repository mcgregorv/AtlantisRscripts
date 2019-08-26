#function for padding out the columns
#used if column headings are say c(1,5,6,10) and want a dataframe with al columns seq(1,10)
pad_cols<-function(data,all_cols,fill=NA){
  
  new_data<-data.frame(matrix(fill,ncol=length(all_cols),nrow=nrow(temp)))
  colnames(new_data)<-all_cols
  this_match<-match(all_cols,colnames(data))
  
  for(i in 1:length(all_cols)){
    this_match<-match(all_cols[i],colnames(data))
    if(!is.na(this_match)){
      new_data[,i]<-data[,this_match]
    }
  }
  return(new_data)
}