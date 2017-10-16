sort_Bin_by_Dec<-function(seq_vector){
  Bin2Dec <- function(x){
    if(x=="Root"){
      dec<-1
    }else{
      x<-paste0("1",x)
      dec<-sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
    }
    dec
  }
  seq_vector[order(unlist(lapply(seq_vector,Bin2Dec)))]
}
