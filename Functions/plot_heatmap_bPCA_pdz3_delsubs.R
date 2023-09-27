plot_heatmap_bPCA_pdz3_delsubs <- function(){
  
  ## isolate data for the domain
  pdz3_bind_df<-scaled_variants_bPCA[scaled_variants_bPCA$domain == "PSD95-PDZ3",]
  
  pdz3_bind_hm_delsub <- ggplot()+
    geom_tile(data=pdz3_bind_df[pdz3_bind_df$type %in% c("delSub_2", "delSub_1", "delSub")], aes(x=Pos,
                                                                                                 y=type,
                                                                                                 fill=scaled_fitness))+
    scale_fill_gradient2(low="red",mid="white",high="blue",midpoint=1,
                         limits=c(min(na.omit(color_scale)),max(na.omit(color_scale))))+
    scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51", "52"))+
    scale_y_discrete(limits=c("delSub_2","delSub_1"),
                     labels=c("delSub_2", "delSub_1"))+
    geom_text(data=pdz3_bind_df[pdz3_bind_df$type %in% c("delSub_2", "delSub_1", "delSub")], aes(x=Pos,
                                                                                                 y=type,
                                                                                                 label=pdz3_bind_df[pdz3_bind_df$type %in% c("delSub_2", "delSub_1", "delSub")]$significant))+
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
  
  return(pdz3_bind_hm_delsub = pdz3_bind_hm_delsub)
}
