plot_m5p_coeff_ins <- function(){
  
ggplot()+
  geom_boxplot(data = top_m5p[top_m5p$coef_name != "(Intercept)",],
               aes(y = s1,
                   x = reorder(as.character(coef_name), as.numeric(s1)),
                   color = "model 5p")) +
  scale_x_discrete(labels=c("structure before: start",
                            "structure after: end",
                            "structure after: c-termini",
                            "structure: loop",
                            "structure: termini",
                            "structure before: n-termini",
                            "position in termini: -4 and <",
                            "structure after: loop",
                            "position in termini: +4 and >",
                            "position in structure: -4 and <",
                            "position in structure: +4 and >",
                            "structure lenght: short", 
                            "position in structure: 3",
                            "position in structure: -3",
                            "structure before: loop",
                            "position in termini: 1",
                            "position in termini: 2",
                            "position in structure: 0",
                            "structure after: AlpaHelix",
                            "structure lenght: medium", 
                            "structure before: AlpaHelix",
                            "position in structure: 2",
                            "position in structure: -2",
                            "position in termini: -1",
                            "position in termini: -3",
                            "position in structure: -1",
                            "position in termini: -2",
                            "structure: Strand",
                            "position in structure: 1",
                            "structure after: Strand",
                            "DDMut",
                            "position in termini: 3",
                            "structure before: Strand",
                            "structure: AlphaHelix"))+
  xlab("predictive features") +
  ylab("coefficients") +
  scale_color_manual(values = c("model 5p" = "yellow3")) +
  theme_classic()+
  theme(axis.text.x = element_text(size=16,
                                   angle = 90, 
                                   hjust = 1),
        axis.text.y = element_text(size=16),
        legend.title = element_blank(),
        legend.text = element_text(size=16),
        axis.title = element_blank())
}
