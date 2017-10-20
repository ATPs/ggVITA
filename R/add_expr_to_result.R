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
  
  #############################################################################
  
  
  
  ggS_ann<-add_expr_2_one_tr(expr_file,
                             fun_alml_readin,
                             result.nb,
                             SorT="S",
                             branch_size=branch_size,
                             expr_size=expr_size,
                             expr_alpha=expr_alpha,
                             tip_size=tip_size,
                             tiplab_size=tiplab_size,
                             mc.cores=mc.cores,
                             colors_gradient=colors_gradient,
                             col_names= col_names)
  
  ggT_ann<-add_expr_2_one_tr(expr_file,
                             fun_alml_readin,
                             result.nb,
                             SorT="T",
                             branch_size=branch_size,
                             expr_size=expr_size,
                             expr_alpha=expr_alpha,
                             tip_size=tip_size,
                             tiplab_size=tiplab_size,
                             mc.cores=mc.cores,
                             colors_gradient=colors_gradient,
                             col_names= col_names)
  
  multiplot(ggS_ann,ggT_ann+scale_x_reverse(),ncol=2)
  
}
  





  
  
 


