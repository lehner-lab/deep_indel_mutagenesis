plot_lenght_fig3 <- function(df,unique_counts){
  
  color1 <- adjustcolor("#414487FF", alpha.f = 0.8) 
  color2 <- adjustcolor("#2A788EFF", alpha.f = 1) 
  color3 <- adjustcolor("#7AD151FF", alpha.f = 0.6)  
  
  min <- min(c(tsuboyama_nat_doms_all$ddG_ML_dels, tsuboyama_nat_doms_all$ddG_ML_ins, tsuboyama_nat_doms_all$ddG_ML_subs))
  max <- max(c(tsuboyama_nat_doms_all$ddG_ML_dels, tsuboyama_nat_doms_all$ddG_ML_ins, tsuboyama_nat_doms_all$ddG_ML_subs))
  
  
  
  ggplot(data = df, 
         aes(x = reorder(as.character(element_lenght), as.numeric(element_lenght)), 
             y = ddG_ML_subs)) +
    geom_boxplot(aes(fill = "ddG_ML_subs"), width = 0.5) +
    geom_jitter(height = 0, width = 0.02) +
    geom_boxplot(data = df, 
                 aes(x = reorder(as.character(element_lenght), as.numeric(element_lenght)), 
                     y = ddG_ML_ins, fill = "ddG_ML_ins"), 
                 width = 0.5, position = position_dodge(width = 0.75)) +
    geom_boxplot(data = df, 
                 aes(x = reorder(as.character(element_lenght), as.numeric(element_lenght)), 
                     y = ddG_ML_dels, fill = "ddG_ML_dels"), 
                 width = 0.5, position = position_dodge(width = 0.75)) +
    scale_fill_manual(values = c("ddG_ML_subs" = color1, "ddG_ML_ins" = color2, "ddG_ML_dels" = color3),
                      labels = c("ddG_ML_subs" = "subs", "ddG_ML_ins" = "ins", "ddG_ML_dels" = "dels")) +  
    scale_y_reverse(limits=c(max,min))+
    geom_text(data = unique_counts, 
              aes(x = reorder(as.character(element_lenght), as.numeric(element_lenght)), 
                  y = -1.9, label = element_no_simple), 
              vjust = 1, size = 5) +
    theme_classic() +
    ylab("ddG") +
    theme(axis.title.x = element_text(size=18),
          axis.text.x = element_text(size=18),
          axis.text.y =  element_text(size=18),
          axis.title.y = element_text(size=18),
          legend.position = "none", 
          legend.title = element_blank())
  
}
