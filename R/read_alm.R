read_alm<-function(file){
  label<-read.table(file,header = T,colClasses = "character") ##################
  label_list<-list()
  label$Lineage<-label$Lineage %>% as.character()
  for(i in 1:nrow(label)){
    label_list[[label[i,1]]]<-label[i,3]
  }
  label_list<<-label_list
  #return(label_list)
  print("The variable label_list is created.")
}