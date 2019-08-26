new.graph<-function (size = 1, width=1, pointsize = 10, xaxs = "i", yaxs = "i", las = 1, quiet=FALSE,filename="", filetype="jpeg", ...) 
{
  windows.options(reset = TRUE)
  width <- width*(21/2.54) * 0.8
  height <- (29.7/2.54) * 0.8
  
  if(quiet==TRUE){
    if(filetype=="svg"){
      svg(filename=filename,
          width = width, height = height*size,  pointsize = pointsize,
          bg = "white", family = "")
    }else{
      jpeg(filename = filename,
           width = width, height = height*size, units = "in", pointsize = pointsize,
           quality = 75,
           bg = "white", res = 700, family = "", restoreConsole = TRUE,
           type = c("windows", "cairo"))
    }
  } else{
    windows(width = width, height = height * size, pointsize = pointsize, 
            bg = "white", rescale = "fixed")
    par(xaxs = xaxs, yaxs = yaxs, las = las, ...)
    .par <<- par(no.readonly = TRUE)
    .par$size <<- size
    .par$pointsize <<- pointsize
    invisible()
  }
 
 
}