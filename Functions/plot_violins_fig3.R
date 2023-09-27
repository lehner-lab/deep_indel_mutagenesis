plot_violins_fig3 <- function(df){
 
  color1 <- adjustcolor("#414487FF", alpha.f = 0.8) 
  color2 <- adjustcolor("#2A788EFF", alpha.f = 0.6) 
  color3 <- adjustcolor("#7AD151FF", alpha.f = 0.6)  
  
  min <- min(c(tsuboyama_nat_doms_all$ddG_ML_dels, tsuboyama_nat_doms_all$ddG_ML_ins, tsuboyama_nat_doms_all$ddG_ML_subs))
  max <- max(c(tsuboyama_nat_doms_all$ddG_ML_dels, tsuboyama_nat_doms_all$ddG_ML_ins, tsuboyama_nat_doms_all$ddG_ML_subs))
  
  
   ggplot(data = df,
         aes(mut, ddG, fill = mut)) +
    geom_violin(position = position_dodge(width = 0.9), alpha = 0.7) +
    geom_jitter(position = position_dodge(width = 0.9), aes(color = mut), alpha = 0.7) +
    scale_fill_manual(values = c(color1, color2, color3)) +
    scale_color_manual(values = c(color1, color2, color3)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size=22),
      axis.text.y = element_text(size=22),
      axis.title = element_blank(),
      legend.title = element_blank(),
      legend.position = "none")+
    scale_y_reverse(limits=c(max,min))
}
