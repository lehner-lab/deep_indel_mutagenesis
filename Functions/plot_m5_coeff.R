plot_m5_coeff <- function(){
  
ggplot()+
  geom_boxplot(data = top_m5[top_m5$coef_name != "(Intercept)",],
               aes(y = s1,
                   x = reorder(as.character(coef_name), as.numeric(s1)),
                   color = "model 5")) +
  scale_x_discrete(labels=c("direct n-neighbour: start", 
                            "structure after: end",
                            "structure: termini",
                            "position in termini: +4 and >",
                            "position in termini: -4 and <",
                            "direct n-neighbour: AlphaHelix",
                            "structure after: ctermini",
                            "structure before: start",
                            "position in structure: 4+",
                            "position in structure: 0",
                            "direct c-neighbour: ctermini",
                            "position in termini: +2",
                            "position in termini: -3",
                            "structure after: AlphaHelix",
                            "direct n-neighbour: Strand",
                            "position in structure: 3",
                            "structure lenght: short", 
                            "position in structure: 2",
                            "position in termini: +1",
                            "position in structure: -4 and <",
                            "structure lenght: medium", 
                            "direct n-neighbour: ntermini",
                            "direct c-neighbour: loop",
                            "position in structure: 1",
                            "structure after: loop",
                            "structure before: loop",
                            "position in structure: -3",
                            "position in structure: -2",
                            "structure after: Strand",
                            "position in termini: -2",
                            "position in termini: -1",
                            "position in structure: -1",
                            "structure: Strand",
                            "direct n-neighbour: loop",
                            "position in termini: +3",
                            "structure before: AlphaHelix",
                            "direct c-neighbour: AlphaHelix",
                            "direct c-neighbour: 310Helix",
                            "direct n-neighbour: 310Helix",
                            "structure: AlphaHelix",
                            "structure before: Strand",
                            "ddG_subs",
                            "direct c-neighbour: Strand"))+
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

