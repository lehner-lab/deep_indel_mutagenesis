plot_m5p_coeff <- function(){
  
ggplot() +
  geom_boxplot(data = top_m5p[top_m5p$coef_name != "(Intercept)",],
               aes(y = s1,
                   x = reorder(as.character(coef_name), as.numeric(s1)),
                   color = "model 5p")) +
  scale_x_discrete(labels=c("direct n-neighbour: start", 
                            "structure after: end",
                            "structure: termini",
                            "position in termini: -4 and <",
                            "position in termini: +4 and >",
                            "direct n-neighbour: AlphaHelix",
                            "position in termini: +2",
                            "structure after: ctermini",
                            "position in termini: -3",
                            "structure before: start",
                            "direct c-neighbour: ctermini",
                            "position in structure: 4+ and >",
                            "position in structure: -4 and <",
                            "position in structure: 0",
                            "position in structure: 2",
                            "direct n-neighbour: Strand",
                            "position in termini: +1",
                            "structure before: ntermini",
                            "structure after: AlphaHelix",
                            "direct c-neighbour: loop",
                            "structure lenght: short",
                            "structure after: loop",
                            "structure lenght: medium", 
                            "position in structure: 3",
                            "structure before: loop",
                            "direct n-neighbour: ntermini",
                            "position in structure: -2",
                            "position in structure: 1",
                            "position in termini: -3",
                            "position in termini: -2",
                            "direct c-neighbour: AlphaHelix",
                            "position in structure: -1",
                            "position in termini: 3",
                            "structure after: Strand",
                            "position in termini: -1",
                            "direct n-neighbour: loop",
                            "structure before: AlphaHelix",
                            "structure: Strand",
                            "direct c-neighbour: 310Helix",
                            "DDMut",
                            "direct n-neighbour: 310Helix",
                            "structure: AlphaHelix",
                            "structure before: Strand",
                            "direct c-neighbour: Strand"))+
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