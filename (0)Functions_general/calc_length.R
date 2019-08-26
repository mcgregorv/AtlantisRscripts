calc_length<-function(w, this_a, this_b){
  ## weight in tonnes
  ((w*1e+6)/this_a)^(1/this_b)
}