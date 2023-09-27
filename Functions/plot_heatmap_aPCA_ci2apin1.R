plot_heatmap_aPCA_ci2apin1 <- function(){
  
  ## isolate data for the domain
  ci2apin1_fold_variants<-scaled_variants_aPCA[scaled_variants_aPCA$domain == "CI2A-PIN1",]
  
  ci2apin1_fold_hm_avg <- ggplot()+
    geom_tile(data=ci2apin1_fold_variants[ci2apin1_fold_variants$type %in% c("delSub_2", "delSub_1", "delSub")], aes(x=Pos, 
                                                                                                         y=type,
                                                                                                         fill=scaled_fitness))+
    geom_tile(data=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "singleDEL",], aes(x=Pos,
                                                                                          y=type,
                                                                                          fill=scaled_fitness))+
    geom_tile(data=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "doubleDEL",], aes(x=Pos, 
                                                                                          y=type,
                                                                                          fill=scaled_fitness))+
    geom_tile(data=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "tripleDEL",], aes(x=Pos, 
                                                                                          y=type,
                                                                                          fill=scaled_fitness))+
    geom_tile(data=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "singleINS",], aes(x=Pos, 
                                                                                          y=type,
                                                                                          fill=scaled_fitness))+
    geom_tile(data=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "doubleINS",], aes(x=Pos, 
                                                                                          y=type,
                                                                                          fill=scaled_fitness))+
    geom_tile(data=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "tripleINS",],aes(x=Pos, 
                                                                                         y=type,
                                                                                         fill=scaled_fitness))+
    scale_fill_gradient2(low="red",mid="white",high="blue",midpoint=1,
                         limits=c(min(na.omit(color_scale)),max(na.omit(color_scale))))+
    scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51", "52","53","54","55","56","57","58","59","60","61","62","63"))+
    scale_y_discrete(limits=c("delSub", "singleDEL", "doubleDEL", "tripleDEL", "singleINS","doubleINS","tripleINS"),
                     labels=c("7","6","5","4","3","2","1"))+
    geom_text(data=ci2apin1_fold_variants[ci2apin1_fold_variants$type %in% c("delSub_2", "delSub_1", "delSub")], aes(x=Pos, 
                                                                                                         y=type,
                                                                                                         label=ci2apin1_fold_variants[ci2apin1_fold_variants$type %in% c("delSub_2", "delSub_1", "delSub")]$significant))+
    geom_text(data=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "singleDEL",], aes(x=Pos, 
                                                                                          y=type,
                                                                                          label=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "singleDEL",]$significant))+
    geom_text(data=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "doubleDEL",], aes(x=Pos, 
                                                                                          y=type,
                                                                                          label=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "doubleDEL",]$significant))+
    geom_text(data=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "tripleDEL",], aes(x=Pos, 
                                                                                          y=type,
                                                                                          label=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "tripleDEL",]$significant))+
    geom_text(data=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "singleINS",], aes(x=Pos, 
                                                                                          y=type,
                                                                                          label=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "singleINS",]$significant))+
    geom_text(data=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "doubleINS",], aes(x=Pos, 
                                                                                          y=type,
                                                                                          label=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "doubleINS",]$significant))+
    geom_text(data=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "tripleINS",], aes(x=Pos, 
                                                                                          y=type,
                                                                                          label=ci2apin1_fold_variants[ci2apin1_fold_variants$type == "tripleINS",]$significant))+
    ylab("CI2A-PIN1")+
    theme_gray()+
    theme(axis.text.x = element_text(size = 13, 
                                     vjust = 0.5, 
                                     hjust = 0.5, 
                                     angle = 0), 
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size=18),
          axis.title.x = element_blank(),
          legend.position = "none")
  
  return(ci2apin1_fold_hm_avg = ci2apin1_fold_hm_avg)
}
