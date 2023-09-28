plot_m5_coeff <- function(){
  
ggplot()+
  geom_boxplot(data = top_m5[top_m5$coef_name != "(Intercept)",],
               aes(y = s1,
                   x = reorder(as.character(coef_name), as.numeric(s1)),
                   color = "model 5")) +
  scale_x_discrete(labels=c("structure after: end", 
                            "structure before: start",
                            "structure: loop",
                            "structure after: c-termini", 
                            "structure: termini",
                            "position in termini: +4 and >",
                            "position in termini: -4 and <",
                            "structure after: loop",
                            "structure before: n-termini",
                            "position in structure: 4+",
                            "position in structure: 0",
                            "position in termini: +2",
                            "position in termini: -3",
                            "position in structure: 3",
                            "position in termini: +1",
                            "structure lenght: short", 
                            "position in structure: 2",
                            "structure before: loop",
                            "position in structure: -4 and <",
                            "structure after: AlphaHelix",
                            "structure lenght: medium", 
                            "position in structure: -1",
                            "position in termini: -1",
                            "position in structure: -3",
                            "position in structure: -2",
                            "position in termini: -2",
                            "structure: Strand",
                            "position in structure: -1",
                            "position in termini: +3",
                            "structure before: AlphaHelix",
                            "structure after: Strand",
                            "structure before: Strand",
                            "structure: AlphaHelix",
                            "ddG_subs"))+
  xlab("predictive features") +
  ylab("coefficients") +
  scale_color_manual(values = c("model 5" = "#C0E8A8", "model 5p" = "yellow3")) +
  theme_classic()+
  theme(axis.text.x = element_text(size=16,
                                   angle = 90, 
                                   hjust = 1),
        axis.text.y = element_text(size=16),
        legend.title = element_blank(),
        legend.text = element_text(size=16),
        axis.title = element_blank())
  
}

