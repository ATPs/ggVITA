Bin2Dec<- function(x){
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}
