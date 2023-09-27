plot_heatmap_bPCA_pdz3_ins <- function(){
  
  ## isolate data for the domain
  pdz3_bind_df<-scaled_variants_bPCA[scaled_variants_bPCA$domain == "PSD95-PDZ3",]
  
  pdz3_bind_hm_avg <- ggplot()+
    geom_tile(data=pdz3_bind_df[pdz3_bind_df$type == "singleINS",], aes(x=Pos, 
                                                                        y=type,
                                                                        fill=scaled_fitness))+
    geom_tile(data=pdz3_bind_df[pdz3_bind_df$type == "doubleINS",], aes(x=Pos, 
                                                                        y=type,
                                                                        fill=scaled_fitness))+
    geom_tile(data=pdz3_bind_df[pdz3_bind_df$type == "doubleINS_CNN",], aes(x=Pos, 
                                                                        y=type,
                                                                        fill=scaled_fitness))+
    geom_tile(data=pdz3_bind_df[pdz3_bind_df$type == "doubleAA",], aes(x=Pos, 
                                                                            y=type,
                                                                            fill=scaled_fitness))+
    geom_tile(data=pdz3_bind_df[pdz3_bind_df$type == "doubleCC",], aes(x=Pos, 
                                                                            y=type,
                                                                            fill=scaled_fitness))+
    geom_tile(data=pdz3_bind_df[pdz3_bind_df$type == "tripleINS",],aes(x=Pos, 
                                                                       y=type,
                                                                       fill=scaled_fitness))+
    geom_tile(data=pdz3_bind_df[pdz3_bind_df$type == "tripleINS_CNNN",],aes(x=Pos, 
                                                                       y=type,
                                                                       fill=scaled_fitness))+
    scale_fill_gradient2(low="red",mid="white",high="blue",midpoint=1,
                         limits=c(min(na.omit(color_scale)),max(na.omit(color_scale))))+
    scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51", "52"))+
    scale_y_discrete(limits=c("singleINS","doubleINS","doubleINS_CNN","doubleAA", "doubleCC", "tripleINS", "tripleINS_CNNN"),
                     labels=c("1x CCC ins", "2x CCC ins", "2x aa repeats","2x Ala ins", "2x Cys ins", "3x CCC ins", "3x aa ins"))+
    geom_text(data=pdz3_bind_df[pdz3_bind_df$type == "singleINS",], aes(x=Pos, 
                                                                        y=type,
                                                                        label=pdz3_bind_df[pdz3_bind_df$type == "singleINS",]$significant))+
    geom_text(data=pdz3_bind_df[pdz3_bind_df$type == "doubleINS",], aes(x=Pos, 
                                                                        y=type,
                                                                        label=pdz3_bind_df[pdz3_bind_df$type == "doubleINS",]$significant))+
    geom_text(data=pdz3_bind_df[pdz3_bind_df$type == "doubleINS_CNN",], aes(x=Pos, 
                                                                        y=type,
                                                                        label=pdz3_bind_df[pdz3_bind_df$type == "doubleINS_CNN",]$significant))+
    geom_text(data=pdz3_bind_df[pdz3_bind_df$type == "doubleAA",], aes(x=Pos, 
                                                                        y=type,
                                                                        label=pdz3_bind_df[pdz3_bind_df$type == "doubleAA",]$significant))+
    geom_text(data=pdz3_bind_df[pdz3_bind_df$type == "doubleCC",], aes(x=Pos, 
                                                                        y=type,
                                                                        label=pdz3_bind_df[pdz3_bind_df$type == "doubleCC",]$significant))+
    geom_text(data=pdz3_bind_df[pdz3_bind_df$type == "tripleINS",], aes(x=Pos, 
                                                                        y=type,
                                                                        label=pdz3_bind_df[pdz3_bind_df$type == "tripleINS",]$significant))+
    geom_text(data=pdz3_bind_df[pdz3_bind_df$type == "tripleINS_CNNN",], aes(x=Pos, 
                                                                        y=type,
                                                                        label=pdz3_bind_df[pdz3_bind_df$type == "tripleINS_CNNN",]$significant))+
    ylab("PSD95-PDZ3")+
    theme_gray()+
    theme(axis.text.x = element_text(size = 13, 
                                     vjust = 0.5, 
                                     hjust = 0.5, 
                                     angle = 0), 
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=14),
          axis.title.x = element_blank(),
          legend.position = "none")
  
  return(pdz3_bind_hm_avg = pdz3_bind_hm_avg)
}
