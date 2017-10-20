add_expr_2_one_tr<-function(expr_file,
                            fun_alml_readin,
                            result.nb,
                            SorT,
                            branch_size=branch_size,
                            expr_size=expr_size,
                            expr_alpha=expr_alpha,
                            tip_size=tip_size,
                            tiplab_size=tiplab_size,
                            mc.cores=mc.cores,
                            colors_gradient=colors_gradient,
                            col_names= col_names
){
  
  
  
  epic_gene_expr<-fread(expr_file)
  
  if(!all(col_names%in% colnames(epic_gene_expr))){
    stop("col_names are not in colnames(the_gene_exprfile)!")
  }
  
  
  epic_gene_expr$Lineage<-epic_gene_expr$cell %>% ggVITA::LN_to_Bin(.)
  
  col_names<-c(col_names,"Lineage")
  
  epic_gene_expr_simple<-data.frame(epic_gene_expr)[,col_names] %>% data.table()
  
  full_tr<-fun_alml_readin$result_list[[result.nb]]  
  
  full_tr2<- ggtree_result(full_tr,
                           isprint = F,
                           branch_size=branch_size,
                           tip_size=tip_size,
                           tiplab_size=tiplab_size)
  
  
  
  
  
  full_tr_nodes_order<-
    full_tr[[paste0("tree",SorT)]]$nodes_order %>% 
    data.table()
  
  setnames(full_tr_nodes_order,"node.order","node")
  
  
  full_tr_ggtree_data<-full_tr2[[paste0("gg",SorT)]]$data %>% data.table()
  
  
  
  full_tr_merge<-
    merge(full_tr_ggtree_data,
          full_tr_nodes_order[,c("parent.seq","parent.order","node.seq","node")],
          by="node"
    )
  
  setkey(full_tr_merge,"node.seq")
  
  
  head(full_tr_merge)
  is.data.table(epic_gene_expr_simple)
  ############################################
  
  
  epic_gene_expr_simple<-epic_gene_expr_simple %>% filter(Lineage %in%  full_tr_merge$node.seq)
  
  
  epic_gene_expr_simple$node.x  <-epic_gene_expr_simple$Lineage %>% mclapply(function(m){full_tr_merge[m]$"x"},mc.cores = mc.cores) %>% unlist()
  
  
  ##
  epic_gene_expr_simple$parent.x<-
    epic_gene_expr_simple$Lineage %>% mclapply(
      function(m){tmp_parent_seq<-as.character(full_tr_merge[m]$parent.seq)
      full_tr_merge[node.seq==tmp_parent_seq,]$"x"
    },
    mc.cores = mc.cores) %>%
    unlist()
  
  ##
  epic_gene_expr_simple<-epic_gene_expr_simple %>% data.table()
  
  ##
  epic_gene_expr_simple_celltime_freq<-
    epic_gene_expr_simple[,"cell"] %>% 
    table() %>% data.table()
  
  
  
  setnames(epic_gene_expr_simple_celltime_freq,c(".","N"),c("cell","time_freq"))
  
  ##
  epic_gene_expr_simple<-
    merge(epic_gene_expr_simple,epic_gene_expr_simple_celltime_freq,by="cell") %>% data.table()
  
  ##
  
  is.data.table(epic_gene_expr_simple)
  
  
  epic_gene_expr_simple$time_rank_in_cell<-epic_gene_expr_simple[,frank(time),by=cell]$V1
  

  print(epic_gene_expr_simple$time_rank_in_cell)
  
  epic_gene_expr_simple<-mutate(epic_gene_expr_simple,seg_x_start=((time_rank_in_cell-1)/ time_freq)*(node.x-parent.x)+parent.x)
  
  
  epic_gene_expr_simple<-mutate(epic_gene_expr_simple,seg_x_end=((time_rank_in_cell)/ time_freq)*(node.x-parent.x)+parent.x)
  
  
  epic_gene_expr_simple$seg_y<-
    epic_gene_expr_simple$Lineage %>% 
    mclapply(function(x){
      full_tr_merge[x]$y
    },mc.cores = mc.cores) %>% 
    unlist()
  
  epic_gene_expr_simple$branch<-
    epic_gene_expr_simple$Lineage %>% 
    mclapply(function(x){
      full_tr_merge[x]$branch
    },mc.cores = mc.cores) %>% 
    unlist()
  
  
  epic_gene_expr_simple$scale_blot<-scale(epic_gene_expr_simple$blot)
  

  
  
  
  
  
  
  
  
  EPIC_colors_gradient<-colors_gradient
  
  ggtr_anotation<-full_tr2[[paste0("gg",SorT)]]+
    geom_segment(aes(x=seg_x_start,
                     xend=seg_x_end,
                     y=seg_y,
                     yend=seg_y,
                     color=scale_blot,
                     group=group
    ),
    linetype="solid",
    size=expr_size,
    alpha=expr_alpha,
    data=epic_gene_expr_simple %>% mutate(group="1"))+
    scale_color_gradientn(colors =EPIC_colors_gradient)+
    geom_tippoint(size=tip_size,aes(fill=I(colorlabel)),shape=21,color="NA")
  
  return(ggtr_anotation)
  
  
}
