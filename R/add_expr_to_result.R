## Date 2017-10-13
## Author : Yuan Meng


## Read the data of expression ,then add to the result

#

add_expr_to_result<-function(expr_file,
                             fun_alml_readin,
                             result.nb,
                             branch_size=0.5,
                             expr_size=0.5,
                             tip_size=0.5,
                             tiplab_size=0.5,
                             mc.cores=16,
                             colors_gradient=c("#CC33FF","#6600FF","#3300FF","#006666","#33FF33","#66FF66","#FFCC00","#FF9900","#FF6600")){
  
  #require("dplyr")
  #require("data.table")
  #require("rlist")
  #require("parallel")
  #require("ggplot2")
  #require("ggtree")
  
  
  
  
  
  add_expr_2_one_tr<-  function(expr_file,fun_alml_readin,result.nb,SorT){
      
      full_tr<-fun_alml_readin$result_list[[result.nb]]
      
      full_tr2<-ggtree_result(full_tr,isprint = F,branch_size=branch_size,tip_size=tip_size,tiplab_size=tiplab_size)
      
      epic_gene_expr<-fread(expr_file)
      
      
      
      if(!all(c("cell","time","blot")%in% colnames(epic_gene_expr))){
        stop("c('cell','time','blot') not in colnames(the_gene_exprfile)!")
      }
      
      epic_gene_expr$Lineage<-epic_gene_expr$cell %>% LN_2_trueLN(.)
      
      
      
      epic_gene_expr_simple<-epic_gene_expr[,c("cell","time","blot","Lineage")]
      
      
      
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
      
      epic_gene_expr_simple$`node.x`<-
        epic_gene_expr_simple$Lineage %>% 
        mclapply(function(x){
          full_tr_merge[x]$`x`
        },mc.cores = mc.cores) %>% 
        unlist()
      
      require(parallel)
      
      epic_gene_expr_simple$`parent.x`<-
        epic_gene_expr_simple$Lineage %>% 
        mclapply(function(x){
          
          tmp_parent_seq<-as.character(full_tr_merge[x]$`parent.seq`)
          
          full_tr_merge[tmp_parent_seq]$`x`
          
        },mc.cores = mc.cores) %>%
        unlist()
      
      epic_gene_expr_simple<-epic_gene_expr_simple %>% data.table()
      
      
      epic_gene_expr_simple_celltime_freq<-
        epic_gene_expr_simple[,"cell"] %>% 
        table() %>% 
        data.table()
      
      
      
      setnames(epic_gene_expr_simple_celltime_freq,c(".","N"),c("cell","time_freq"))
      
      
      epic_gene_expr_simple<-
        merge(epic_gene_expr_simple,epic_gene_expr_simple_celltime_freq,by="cell")
      
      
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
      
      #epic_gene_expr_simple$`blot`[epic_gene_expr_simple$`blot`<0]<-0
      
      
      epic_gene_expr_simple$`scale_blot`<-with(epic_gene_expr_simple,scale(blot))
      
      
      
      epic_gene_expr_simple_na_rm<-
        data.frame(node=epic_gene_expr_simple$node.x,epic_gene_expr_simple) %>% 
        filter(Lineage %in%  full_tr_merge$node.seq)
      
      EPIC_colors_gradient<-colors_gradient
    
      
      #EPIC_colors_gradient<-rainbow(4)
      
      ggtr_anotaiton<-full_tr2[[paste0("gg",SorT)]]+
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
      
      return(ggtr_anotaiton)
      
      
    }
  
  
#############################################################################
  
  
  
    ggS_ann<-add_expr_2_one_tr(expr_file,fun_alml_readin,result.nb,"S")

    ggT_ann<-add_expr_2_one_tr(expr_file,fun_alml_readin,result.nb,"T")
  
  return(multiplot(ggS_ann,ggT_ann+scale_x_reverse(),ncol=2))
  
}
  
  

  
  
  
 

