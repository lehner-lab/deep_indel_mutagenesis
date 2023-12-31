plot_heatmap_bPCA_pdz3_allins <- function(){
  
  pdz3_bind_df <- scaled_variants_bPCA[scaled_variants_bPCA$domain == "PSD95-PDZ3",]
  
  pdz3_bind_allins_hm<-ggplot()+
    geom_tile(data=pdz3_bind_df[pdz3_bind_df$type == "allins",], 
              aes(x=Pos, 
                  y=insID,
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
    scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52"))+
    ylab("pdz3_bind")+
    geom_text(data=pdz3_bind_df[pdz3_bind_df$type == "allins",], 
              aes(x=Pos,
                  y=insID,
                  label=pdz3_bind_df[pdz3_bind_df$type == "allins",]$significant))+
    theme(axis.title.y = element_text(size=18),
          axis.text.y = element_text(size=16))+
    theme_gray()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
  
  return(list(pdz3_bind_allins_hm = pdz3_bind_allins_hm))
  
}
