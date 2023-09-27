plot_heatmap_aPCA_grb2_ins <- function(){
  
  grb2_fold_df <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "GRB2-SH3",]
  
  grb2_fold_hm_avg<-ggplot()+
    geom_tile(data=grb2_fold_df[grb2_fold_df$type == "singleINS",], aes(x=Pos, 
                                     y=type,
                                     fill=scaled_fitness))+
    geom_tile(data=grb2_fold_df[grb2_fold_df$type == "doubleINS",], aes(x=Pos, 
                                      y=type,
                                      fill=scaled_fitness))+
    geom_tile(data=grb2_fold_df[grb2_fold_df$type == "doubleINS_CNN",], aes(x=Pos, 
                                      y=type,
                                      fill=scaled_fitness))+
    geom_tile(data=grb2_fold_df[grb2_fold_df$type == "doubleAA",], aes(x=Pos, 
                                                                    y=type,
                                                                    fill=scaled_fitness))+
    geom_tile(data=grb2_fold_df[grb2_fold_df$type == "doubleCC",], aes(x=Pos, 
                                                                    y=type,
                                                                    fill=scaled_fitness))+
    geom_tile(data=grb2_fold_df[grb2_fold_df$type == "tripleINS",],aes(x=Pos, 
                                      y=type,
                                      fill=scaled_fitness))+
    geom_tile(data=grb2_fold_df[grb2_fold_df$type == "tripleINS_CNNN",], aes(x=Pos, 
                                       y=type,
                                       fill=scaled_fitness))+
    scale_fill_gradient2(low="red",mid="white",high="blue",midpoint=1,
                         limits=c(min(na.omit(color_scale)),max(na.omit(color_scale))))+
    theme(axis.text.x = element_text(size = 11, 
                                     vjust = 0.5, 
                                     hjust = 0.5, 
                                     angle = 0), 
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51", "52"))+
    scale_y_discrete(limits=c("singleINS","doubleINS","doubleINS_CNN", "doubleAA", "doubleCC","tripleINS","tripleINS_CNNN"),
                     labels=c("1x aa repeat", "2x CCC ins", "2x aa repeat", "2x insAA", "2x insCC", "3x CCC ins", "3x aa repeat"))+
    geom_text(data=grb2_fold_df[grb2_fold_df$type == "singleINS",], aes(x=Pos, 
                                     y=type,
                                     label=grb2_fold_df[grb2_fold_df$type == "singleINS",]$significant))+
    geom_text(data=grb2_fold_df[grb2_fold_df$type == "doubleINS",], aes(x=Pos, 
                                      y=type,
                                      label=grb2_fold_df[grb2_fold_df$type == "doubleINS",]$significant))+
    geom_text(data=grb2_fold_df[grb2_fold_df$type == "doubleINS_CNN",], aes(x=Pos, 
                                      y=type,
                                      label=grb2_fold_df[grb2_fold_df$type == "doubleINS_CNN",]$significant))+
    geom_text(data=grb2_fold_df[grb2_fold_df$type == "doubleAA",], aes(x=Pos, 
                                                                    y=type,
                                                                    label=grb2_fold_df[grb2_fold_df$type == "doubleAA",]$significant))+
    geom_text(data=grb2_fold_df[grb2_fold_df$type == "doubleCC",], aes(x=Pos, 
                                                                    y=type,
                                                                    label=grb2_fold_df[grb2_fold_df$type == "doubleCC",]$significant))+
    geom_text(data=grb2_fold_df[grb2_fold_df$type == "tripleINS",], aes(x=Pos, 
                                       y=type,
                                       label=grb2_fold_df[grb2_fold_df$type == "tripleINS",]$significant))+
    geom_text(data=grb2_fold_df[grb2_fold_df$type == "tripleINS_CNNN",], aes(x=Pos, 
                                       y=type,
                                       label=grb2_fold_df[grb2_fold_df$type == "tripleINS_CNNN",]$significant))+
    theme_gray()+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size=14),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.title.x = element_blank(),
          legend.position = "none")
  
  return(list(grb2_fold_hm_avg = grb2_fold_hm_avg))
}
  