plot_heatmap_aPCA_pdz3_subs <- function(){
  
  # isolate PSD95-PDZ3 data
  pdz3_fold_df <- scaled_variants_aPCA[scaled_variants_aPCA$domain == "PSD95-PDZ3",]
  
  pdz3_fold_subs_hm<-ggplot()+
    geom_tile(data=pdz3_fold_df[pdz3_fold_df$mut_type == "substitutions"], 
              aes(x=Pos, 
                  y=Mut,
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
    ylab("pdz3_fold")+
    geom_text(data=pdz3_fold_df[pdz3_fold_df$mut_type == "substitutions"], 
              aes(x=Pos,
                  y=Mut,
                  label=pdz3_fold_df[pdz3_fold_df$mut_type == "substitutions"]$significant))+
    theme(axis.title.y = element_text(size=18),
          axis.text.y = element_text(size=20))+
    theme_gray()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  
  ## make a color code that colors the sequence by secondary struc. 
  pdz3_seq_code <- pdz3_fold_df[,c("Pos","Structure_name")]
  pdz3_seq_code <- pdz3_seq_code[!duplicated(pdz3_seq_code),]
  pdz3_seq_code <- pdz3_seq_code[-1,]
  
  pdz3_seq_code <- pdz3_seq_code %>% mutate(color_seq =
                                              ifelse(Structure_name == "Coil" | Structure_name == "Turn" | Structure_name == "Bridge", "black",
                                                     ifelse(Structure_name == "AlphaHelix" | Structure_name == "310Helix", "red",
                                                            ifelse(Structure_name == "Strand", "blue", NA))))
  
  ## add the last row as strand
  pdz3_seq_code[nrow(pdz3_seq_code),]$color_seq <- "blue"
  
  
  return(list(pdz3_fold_subs_hm = pdz3_fold_subs_hm,
              pdz3_seq_code = pdz3_seq_code))  
}



