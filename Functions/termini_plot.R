termini_plot <- function(df, ddG, unique_counts){
  
  min=min(c(ntermini_ddG$ddG_ML_dels,
            ntermini_ddG$ddG_ML_ins,
            ntermini_ddG$ddG_ML_subs))
  
  max=max(c(ntermini_ddG$ddG_ML_dels,
            ntermini_ddG$ddG_ML_ins,
            ntermini_ddG$ddG_ML_subs))
  
  short<-ggplot(data = df, 
                aes(x = reorder(as.character(align_to_center), as.numeric(align_to_center)), 
                    y = ddG)) +
    geom_boxplot() +
    geom_jitter(height = 0, width = 0.02) +
    scale_y_reverse(limits = c(max, min)) +
    geom_text(data = unique_counts, 
              aes(x = reorder(as.character(align_to_center), as.numeric(align_to_center)), 
                  y = -1.3, label = element_no_simple), 
              vjust = 1, size = 4)+
    theme_classic()+
    ylab("ddG deletions")+
    theme(axis.title.y = element_blank(),
          axis.text.x = element_text(size=14),
          axis.text.y=element_blank(),
          axis.title.x = element_blank())
  
  return(list(short))
  
}
