## Date 2017-10-13
## Author : Yuan Meng


## Read the data of expression ,then add to the result

#

add_expr_to_result<-function(expr_file,
                             fun_alml_readin,
                             result.nb,
                             branch_size=0.5,
                             expr_size=0.5,
                             expr_alpha=0.5,
                             tip_size=0.5,
                             tiplab_size=0.5,
                             mc.cores=16,
                             colors_gradient=c("#CC33FF",
                                               "#6600FF",
                                               "#3300FF",
                                               "#006666",
                                               "#33FF33",
                                               "#66FF66",
                                               "#FFCC00",
                                               "#FF9900",
                                               "#FF6600"),
                             col_names=c("cell","time","blot")){
  
  #require("dplyr")
  #require("data.table")
  #require("rlist")
  #require("parallel")
  #require("ggplot2")
  #require("ggtree")
  
  epic_gene_expr<-data.table::fread(expr_file)
  
  
  if(!all(col_names%in% colnames(epic_gene_expr))){
    stop("col_names are not in colnames(the_gene_exprfile)!")
  }
  
  
  epic_gene_expr$Lineage<-epic_gene_expr$cell %>% LN_to_Bin(.)
  col_names<-c(col_names,"Lineage")
  epic_gene_expr_simple<-data.frame(epic_gene_expr)[,col_names] %>% data.table()
  
  
  
  
  
  
  
  
  
  add_expr_2_one_tr<-function(fun_alml_readin,
                              result.nb,
                              SorT){
    
    full_tr<-fun_alml_readin$result_list[[result.nb]]
    
    
    
    #   think about the args passing to ggtree_result
    environment(ggtree_result)=environment()
    
    ggtree_result2<-ggtree_result
    full_tr2<- ggtree_result2(full_tr,
                                     isprint = F,
                                     branch_size=branch_size,
                                     tip_size=tip_size,
                                     tiplab_size=tiplab_size)
    
    
    full_tr_nodes_order<-
      full_tr[[paste0("tree",SorT)]]$nodes_order %>% 
      data.table()
    
    setnames(full_tr_nodes_order,"node.order","node")
    
    
    full_tr_ggtree_data<-
      full_tr2[[paste0("gg",SorT)]]$data %>% 
      data.table()
    
   
    
    full_tr_merge<-
      merge(full_tr_ggtree_data,
            full_tr_nodes_order[,c("parent.seq","parent.order","node.seq","node")],
            by="node"
      )
    
    setkey(full_tr_merge,"node.seq")
    
    epic_gene_expr_simple$"node.x"  <-
      epic_gene_expr_simple$Lineage %>% 
      mclapply( function(m){full_tr_merge[m]$"x"} , mc.cores = mc.cores ) %>% 
      unlist()
    
    
    ##
    epic_gene_expr_simple$"parent.x"<-
      epic_gene_expr_simple$"Lineage" %>%
      mclapply(function(m){tmp_parent_seq<-as.character(full_tr_merge[m]$`parent.seq`);
      full_tr_merge[tmp_parent_seq]$"x"
      },
      mc.cores = mc.cores) %>%
      unlist()
    
    ##
    epic_gene_expr_simple<-epic_gene_expr_simple %>% data.table()
    
    ##
    epic_gene_expr_simple_celltime_freq<-
      epic_gene_expr_simple[,"cell"] %>% 
      table() %>% 
      data.table()
    
    
    
    setnames(epic_gene_expr_simple_celltime_freq,c(".","N"),c("cell","time_freq"))
    
    ##
    epic_gene_expr_simple<-
      merge(epic_gene_expr_simple,epic_gene_expr_simple_celltime_freq,by="cell")
    
    ##
    epic_gene_expr_simple$`time_rank_in_cell`<- 
      epic_gene_expr_simple[,rank(time),by=cell]$`V1`
    
    
    epic_gene_expr_simple$`seg_x_start`<-with(epic_gene_expr_simple,((time_rank_in_cell-1)/ time_freq)*(node.x-parent.x)+parent.x)
    
    
    epic_gene_expr_simple$`seg_x_end`<-with(epic_gene_expr_simple,((time_rank_in_cell)/ time_freq)*(node.x-parent.x)+parent.x)
    
    
    epic_gene_expr_simple$`seg_y`<-
      epic_gene_expr_simple$Lineage %>% 
      mclapply(function(x){
        full_tr_merge[x]$`y`
      },mc.cores = mc.cores) %>% 
      unlist()
    
    epic_gene_expr_simple$`branch`<-
      epic_gene_expr_simple$Lineage %>% 
      mclapply(function(x){
        full_tr_merge[x]$`branch`
      },mc.cores = mc.cores) %>% 
      unlist()
    
    
    epic_gene_expr_simple$scale_blot<-with(epic_gene_expr_simple,scale(blot))
    
   
    
    
    
    
    
    epic_gene_expr_simple_na_rm<-
      data.frame(node=epic_gene_expr_simple$node.x,epic_gene_expr_simple) %>% 
      filter(Lineage %in%  full_tr_merge$node.seq)
    
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
      data=epic_gene_expr_simple_na_rm %>% mutate(group="1"))+
      scale_color_gradientn(colors =EPIC_colors_gradient)+
      geom_tippoint(size=tip_size,aes(fill=I(colorlabel)),shape=21,color="NA")
    
    return(ggtr_anotation)
    
    
  }
  
  
  #############################################################################
  
  
  
  ggS_ann<-add_expr_2_one_tr(fun_alml_readin,result.nb,SorT="S")
  
  ggT_ann<-add_expr_2_one_tr(fun_alml_readin,result.nb,SorT="T")
  
  multiplot(ggS_ann,ggT_ann+scale_x_reverse(),ncol=2)
  
}
  

  
  
  
 


