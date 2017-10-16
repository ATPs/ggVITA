read_alml<-function(file=""){
  the_result<-list()
  the_prefix<-c("Score","RootS","RootT","PruneS","PruneT","MatchS","MatchT","}","PValue","Min")
  the_text<-as.list(readLines(file))
  nl<-length(the_text)
  for(i in 1:nl){
    the_line<-the_text[[i]]
    if(regexpr("^[0-9]",the_line)){
      num2<-as.integer(regmatches(the_line,regexpr("^([0-9]+)",the_line)))
      if(length(num2)>0){
        num<-as.character(regmatches(the_line,regexpr("^([0-9]+)",the_line)))
        the_result[[num]]<-list()
        the_result[[num]]<-list("num"=c(num))
      }
    }
    for(i2 in 1:(length(the_prefix)-3)){
      if(startsWith(the_line,the_prefix[i2])){
        the_result[[num]][[as.character(the_prefix[i2])]]<-unlist(strsplit(the_line,split = ":"))[2]
      } 
    }
    if(startsWith(the_line,"PValue")){
      the_result[["PValue"]][["all"]]<-unlist(strsplit(the_line,split = ":"))[2]
    }
    if(startsWith(the_line,"Min")){
      the_result[["PValue"]][["Min"]]<-unlist(strsplit(unlist(strsplit(the_line,split = " "))[1],split = ":"))[2]
      the_result[["PValue"]][["Max"]]<-unlist(strsplit(unlist(strsplit(the_line,split = " "))[2],split = ":"))[2]
      the_result[["PValue"]][["AVG"]]<-unlist(strsplit(unlist(strsplit(the_line,split = " "))[3],split = ":"))[2]
    }
  }
  
  result2phylo <- function(the_result) {
    
    # transform a result(S and T) into 2 phylos
    
    # @some functions to use
    
    get_parent <- function(x) {
      x <- substr(x, 1, (nchar(x) - 1))
      if (x == "") {
        x <- "Root"
      }
      x
    }
    
    nb.branches <- function(the_seq, x) {
      
      if(x=="Root"){
        x<-""
      }
      
      x_0 <- paste0(x, "0")
      x_1 <- paste0(x, "1")
      i <- 0
      
      if (x_0 %in% the_seq) {
        i <- i + 1
      }
      if (x_1 %in% the_seq) {
        i <- i + 1
      }
      
      i
    }
    
    
    get_outside_son <- function(the_seq) {
      unique(unlist(lapply(the_seq, function(x) {
        if (nb.branches(the_seq, x) == 1) {
          x_0 <- paste0(x, "0")
          x_1 <- paste0(x, "1")
          return(setdiff(c(x_0, x_1), the_seq))
        }
      })))
    }
    
    get_outside_parent <- function(the_seq) {
      the_seq<-the_seq[-1]
      min_level <- min(nchar(the_seq))
      y <- unlist(lapply(the_seq, function(x) {
        
        if (nchar(x) > min_level & x != "Root" ) {
          if (T) {
            if (!(substr(x, 1, (nchar(x) - 1)) %in% the_seq)) {
              return(substr(x, 1, (nchar(x) - 1)))
            }
          }
        }
      }))
      unique(y)
    }
    
    
    s2v <- function(x) {
      unlist(strsplit(x, split = " "))
    }
    
    
    
    
    
    tr2phylo<-function(seq,the_root){
      
      # @ group the tips and the nodes
      
      the_seq <- s2v(seq)
      
      the_ori_seq_nb <- length(the_seq)
      
      #sub_seq<-unlist(lapply(seq,get_sub_seq))
      
      
      # @ test
      
      if (the_root != "Root") {
        if (!all(startsWith(the_seq, prefix = as.character(the_root)))) {
          stop("Existing nodes outside the subroot!")
        }
      }
      
      
      
      # @ make-up the tree
      
      i = 1
      
      #the_seq<-s2v(fun.alml[[3]]$MatchS)
      
      the_outside_son <- get_outside_son(the_seq)
      
      the_outside_parent <- get_outside_parent(the_seq)
      
      the_outside_node <- union(the_outside_parent, the_outside_son)
      
      
      
      #print(paste0("outside parent: ", paste(the_outside_parent, collapse = " ")))
      
      #print(paste0("outside son: ", paste(the_outside_son, collapse = " ")))
      
      # print(paste0("son==parent: ", paste(
      #   intersect(the_outside_parent, the_outside_son), collapse = " "
      # )))
      
      
      ##################
      
      the_seq <- c(the_seq, the_outside_node)
      
      ##################
      
      the_outside_son <- get_outside_son(the_seq)
      
      the_outside_parent <- get_outside_parent(the_seq)
      
      the_outside_node <- union(the_outside_parent, the_outside_son)
      
      
      # 
      # print(paste0("outside parent: ", paste(the_outside_parent, collapse = " ")))
      # 
      # print(paste0("outside son: ", paste(the_outside_son, collapse = " ")))
      # 
      # print(paste0("son==parent: ", paste(
      #   intersect(the_outside_parent, the_outside_son), collapse = " "
      # )))
      
      
      
      
      
      while (length(the_outside_parent) > 0 |
             length(the_outside_son) > 0) {
        i = i + 1
        
        #print(paste0("make-up times :", i))
        
        the_outside_son <- get_outside_son(the_seq)
        
        the_outside_parent <- get_outside_parent(the_seq)
        
        the_outside_node <- union(the_outside_parent, the_outside_son)
        
        
        # print(paste0("the outside nodes:", paste(the_outside_node, collapse = " ")))
        # 
        # print(paste0("outside parent: ", paste(the_outside_parent, collapse = " ")))
        # 
        # print(paste0("outside son: ", paste(the_outside_son, collapse = " ")))
        # 
        # print(paste0("son==parent: ", paste(
        #   intersect(the_outside_parent, the_outside_son),
        #   collapse = " "
        # )))
        # 
        ##################
        
        the_seq <- c(the_seq, the_outside_node)
        
        ##################
        
        the_outside_son <- get_outside_son(the_seq)
        
        the_outside_parent <- get_outside_parent(the_seq)
        
        the_outside_node <- union(the_outside_parent, the_outside_son)
        
        
        
        # print(paste0("outside parent: ", paste(the_outside_parent, collapse = " ")))
        # 
        # print(paste0("outside son: ", paste(the_outside_son, collapse = " ")))
        # 
        # print(paste0("son==parent: ", paste(
        #   intersect(the_outside_parent, the_outside_son),
        #   collapse = " "
        # )))
        
      }
      
      

      
      the_parent <-unlist(lapply(the_seq[-1],get_parent))
      the_parent <-c(the_seq[1],the_parent)
      
     
      if (any(unlist(lapply(the_seq, function(x) {
        nb.branches(the_seq, x)
      })) == 1)) {
        the_single_banche_node <-
          the_seq[which(unlist(lapply(the_seq, function(x) {
            nb.branches(the_seq, x)
          })) == 1)]
        print(paste0(
          "existing nb.branches:",
          paste(the_single_banche_node, collapse = " ")
        ))
        stop("existing nb.branches==1")
      }
      
      
      
      the_tip <-
        the_seq[which(unlist(lapply(the_seq, function(x) {
          nb.branches(the_seq, x)
        })) == 0)]
      
      the_node<-
        the_seq[which(unlist(lapply(the_seq, function(x) {
          nb.branches(the_seq, x)
        })) == 2)]
      
      
      # @ give the nb.node to the seq
      
      
      
      get_the_order <- function(x) {
        if (x %in% the_tip) {
          x.order <- which(the_tip == x)
        }
        if (x %in% the_node) {
          x.order <- which(the_node == x) + length(the_tip)
        }
        x.order
      }
      
      # @ give the nb.node to the seq
      
      
      
      the_seq_order <- unlist(lapply(the_seq, get_the_order))
      
      the_parent_order <- unlist(lapply(the_parent, get_the_order))
      
      
      nodes_order <- data.frame(
        parent.seq = the_parent,
        parent.order = the_parent_order,
        node.seq = the_seq,
        node.order = the_seq_order
      ,stringsAsFactors = F)
      
      
      
      
      get_the_label <- function(x) {
        x<-as.character(x)
        if (x %in% names(label_list)) {
          x.label <- label_list[[x]]
        }else{
          x.label <- "???"
        }
        x.label
        
      }
      nodes_order<-nodes_order %>% mutate(label=
                                            unlist(lapply(nodes_order$node.seq,function(x){
                                              x<-as.character(x)
                                              if (x %in% names(label_list)) {
                                                x.label <- label_list[[x]]
                                                if (length(x.label) == 0) {
                                                  x.label <- "???"
                                                }
                                              } else{
                                                x.label <- "???"
                                              }
                                              x.label
                                            })))
      
      nodes_order<-nodes_order %>% mutate(isTip=
                                            unlist(lapply(nodes_order$node.seq,function(x){
                                              x<-as.character(x)
                                              if(nb.branches(the_seq,x)==2){
                                                y<-"FALSE"
                                              }
                                              if(nb.branches(the_seq,x)==0){
                                                y<-"TRUE"
                                              }
                                              if(nb.branches(the_seq,x)==1){
                                                print("Error:nb.branches==1")
                                              }
                                              y
                                            })))
      nodes_order <-
        data.table(nodes_order[order(nodes_order$node.order), ])
      
      # @ phylo
      
      phylo_tree <- list()
      
      phylo_tree$edge <-
        matrix(cbind(the_parent_order, the_seq_order), ncol = 2)[-1,]
      
      phylo_tree$tip.label <-
        as.character(unlist(lapply(the_seq[order(the_seq_order)][1:length(the_tip)], get_the_label)))
      
      phylo_tree$Nnode <-
        as.integer(length(unique(the_node)))
      
      class(phylo_tree) <- "phylo"
      
      return(list(phylo_tr = phylo_tree, nodes_order = nodes_order))
      
      
    }## end tr2phylo function
    
    
    #print("===Transforming S to phylo===")
    
    if(is.na(the_result$PruneS)==F){
      treeS<-tr2phylo(paste0(the_result$RootS," ",the_result$MatchS,the_result$PruneS),the_result$RootS)
      
    }else{
      treeS<-tr2phylo(paste0(the_result$RootS," ",the_result$MatchS),the_result$RootS)
    }
    
    treeS$nodes_order<-treeS$nodes_order %>% mutate(matched_or_pruned=
                                                      unlist(lapply(treeS$nodes_order$node.seq,function(m){
                                                        m<-as.character(m)
                                                        if(m %in% s2v(paste0(the_result$RootS," ",the_result$MatchS))){
                                                          p_or_s<-"matched"
                                                        }else if(m %in% s2v(the_result$PruneS)){
                                                          p_or_s<-"pruned"
                                                        }else{
                                                          p_or_s<-"pruned_sister"
                                                        }
                                                        p_or_s
                                                      })))
    
    #print("===Transforming T to phylo===")
    if(is.na(the_result$PruneT)==F){
      treeT<-tr2phylo(paste0(the_result$RootT," ",the_result$MatchT,the_result$PruneT),the_result$RootT)
      
    }else{
      treeT<-tr2phylo(paste0(the_result$RootT," ",the_result$MatchT),the_result$RootT)
    }
    
    treeT$nodes_order<- treeT$nodes_order%>% mutate(matched_or_pruned=
                                                      unlist(lapply(treeT$nodes_order$node.seq,function(m){
                                                        m<-as.character(m)
                                                        if(m %in% s2v(paste0(the_result$RootT," ",the_result$MatchT))){
                                                          p_or_s<-"matched"
                                                        }else if(m %in% s2v(the_result$PruneT)){
                                                          p_or_s<-"pruned"
                                                        }else{
                                                          p_or_s<-"pruned_sister"
                                                        }
                                                        p_or_s
                                                      })))
    
    matched_pair<-data.table(mtS=s2v(paste0(the_result$RootS," ",the_result$MatchS)),mtT=s2v(paste0(the_result$RootT," ",the_result$MatchT)))
    
    
    return(list(treeS=treeS,treeT=treeT,matched_pair=matched_pair,the_result=the_result))
    
  }
  
  result_list<-lapply(the_result[1:(length(the_result)-1)],function(x){print(x$num);return(result2phylo(x))})
  final_result<-list(result_list=result_list,result_analysis=the_result[length(the_result)])
  return(final_result)
}
