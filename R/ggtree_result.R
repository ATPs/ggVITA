ggtree_result<-function(one_result,layout="rectangular",branch_size=0.2,tip_size=0.5,tiplab_size=0.1,isprint=TRUE){
 
  
  
  #require("dplyr")
  #require("data.table")
  #require("rlist")
  #require("parallel")
  #require("ggplot2")
  #require("ggtree")
  #require("pipeR")
  
 
  s2v<-function(x){
    unlist(strsplit(x,split=" "))
  }
  
  
  S_pruned<-filter(one_result$treeS$nodes_order,matched_or_pruned=="pruned")$node.order
  T_pruned<-filter(one_result$treeT$nodes_order,matched_or_pruned=="pruned")$node.order
 
  
  
  trS<-one_result$treeS$phylo_tr %>% groupOTU(focus=c(1)) %>% ggtree(aes(linetype=group),size=branch_size,layout = layout,alpha=1)
  trS$data$group[T]<-0

  if(length(S_pruned)!=0){
    trS$data$group[trS$data$node %in% S_pruned]<-1
  }
  
 
  trT<-one_result$treeT$phylo_tr %>% groupOTU(focus=c(1)) %>% ggtree(aes(linetype=group),size=branch_size,layout = layout,alpha=1)
  trT$data$group[T]<-0
  
  if(length(T_pruned)!=0){
    trT$data$group[trT$data$node %in% T_pruned]<-1
  }
  
  
  ### create the ggtree data
  dt_S<-merge(trS$data,one_result$treeS$nodes_order[,c("parent.seq","parent.order","node.seq","node.order","matched_or_pruned")],by.x="node",by.y="node.order")
  #dt_S %>% View()
  dt_T<-merge(trT$data,one_result$treeT$nodes_order[,c("parent.seq","parent.order","node.seq","node.order","matched_or_pruned")],by.x="node",by.y="node.order")
  #dt_T %>% View()
  
  matched_S_and_T<-list.parse(one_result$matched_pair)
  
  
  # pruned_sister
  for(i in list.parse(filter(dt_S,matched_or_pruned=="pruned_sister"))){
    if(T){
      dt_S[dt_S$node.seq==i$parent.seq,"y"]<-dt_S[dt_S$node.seq==i$node.seq,"y"]
    }
  }
  
  for(i in matched_S_and_T){
    dt_T[dt_T$node.seq==i$mtT,"y"]<-dt_S[dt_S$node.seq==i$mtS,"y"]
  }
  
  for(i in list.parse(filter(dt_T,matched_or_pruned=="pruned_sister"))){
    dt_T[dt_T$node.seq==i$node.seq,"y"]<-dt_T[dt_T$node.seq==i$parent.seq,"y"]
  }
  
  for(i in list.parse(filter(dt_T,matched_or_pruned=="pruned"))%>>%list.sort(parent.order)){
    sub_tips_seq<-dt_T$node.seq[unlist(lapply(dt_T$node.seq,
                                      function(x){
                                        startsWith(x,i$parent.seq)
                                        }
                                      )
                               )
                        ]
    
    sub_tips_seq<-setdiff(sub_tips_seq,i$node.seq)
    subtips_y_floor<-min(filter(dt_T,node.seq %in% sub_tips_seq)$y)
    subtips_y_ceiling<-max(filter(dt_T,node.seq %in% sub_tips_seq)$y)
    
    sister.seq<-paste0(substr(i$node.seq,1,nchar(i$node.seq)-1),
                       list("1"="0","0"="1")[[substr(i$node.seq,nchar(i$node.seq),nchar(i$node.seq))]])
    
    
    if(i$isTip==T & filter(dt_T,node.seq==sister.seq)$isTip==T){
      if(substr(i$node.seq,nchar(i$node.seq),nchar(i$node.seq))=="0"){
        dt_T[dt_T$node.seq==i$node.seq,"y"]<- dt_T[dt_T$node.seq==i$parent.seq,"y"]-0.35
      }else{
        dt_T[dt_T$node.seq==i$node.seq,"y"]<- dt_T[dt_T$node.seq==i$parent.seq,"y"]+0.35
      }
    }else{
      if(substr(i$node.seq,nchar(i$node.seq),nchar(i$node.seq))=="0"){
        dt_T[dt_T$node.seq==i$node.seq,"y"]<-subtips_y_floor-0.35
      }else{
        dt_T[dt_T$node.seq==i$node.seq,"y"]<-subtips_y_ceiling+0.35
      }
    }
}
  
  Bin2Dec <- function(x) {
    if (x == "Root") {
      dec <- 1
    } else{
      x <- paste0("1", x)
      dec <-
        sum(2 ^ (which(rev(
          unlist(strsplit(as.character(x), "")) == 1
        )) - 1))
    }
    dec
  }
  
  
  Dec2Bin <- function(x) {
    if (x == 1) {
      bin <- "Root"
    } else{
      i <- 0
      string <- numeric(32)
      while (x > 0) {
        string[32 - i] <- x %% 2
        x <- x %/% 2
        i <- i + 1
      }
      first <- match(1, string)
      bin <-
        sub("1", "", paste(as.character(string[first:32]), collapse = ""))
    }
    bin
  }
  
  
  
  #tmp[which((tmp-99)==min(abs(tmp-99)))] %>% Dec2Bin(.)
  
  
  
  
  
  # dt_T$branch<-dt_T$branch/1.5
  # dt_S$branch<-dt_S$branch/1.5
  
  
  trT$data<- merge(trT$data[,-4],dt_T[,c("node","y")],by="node")
  trS$data<- merge(trS$data[,-4],dt_S[,c("node","y")],by="node")
  
  
  
  ###############################
  ## put two tree together and color the tips by tissue classes
  
  labellist2<-
    c(
    "Neu" = c("blue"),
    "Dea" = c("green"),
    "Str" = c("red"),
    "Epi" = c("orange"),
    "Mus" = c("yellow"),
    "Bla" = c("purple"),
    "Gla" = c("black"),
    "Int" = c("brown"),
    "Ger" = c("cyan"),
    "???" = c("skyblue")
  ) 
  
  
  trS$data$colorlabel<-labellist2[trS$data$label] 
  trT$data$colorlabel<-labellist2[trT$data$label]
  
  ggS<-trS+
    geom_tippoint(size=tip_size,aes(fill=I(colorlabel)),shape=21,color="NA")+
    geom_tiplab(align = T,size=tiplab_size)+
    ggtitle(paste0("result.nb:",one_result$the_result$num,"  ","score=",one_result$the_result$Score))
  ggT<-trT+
    scale_x_reverse()+
    geom_tippoint(size=tip_szie,aes(fill=I(colorlabel)),shape=21,color="NA")+
    geom_tiplab(align = T,size=tiplab_size,hjust =1 )+
    ggtitle(paste0("RootS:",one_result$the_result$RootS,"  ","RootT:",one_result$the_result$RootT))
  
  if(isprint==T){multiplot(ggS,ggT,ncol = 2)}
  
  return(list("ggS"=ggS,"ggT"=ggT))

}
