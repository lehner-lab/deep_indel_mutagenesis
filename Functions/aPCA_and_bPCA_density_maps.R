aPCA_and_bPCA_density_maps <- function(df_bind, df_fold){

ggplot() +
  geom_density_ridges(data = df_bind,
                      aes(x = scaled_fitness, y = mut_type, color = mut_type, fill = mut_type)) +
  geom_density_ridges(data = df_fold,
                      aes(x = scaled_fitness, y = mut_type, color = "abundance", fill = "abundance"),
                      scale = 1.6,
                      alpha = 0.3) +
  geom_vline(aes(xintercept = 1),
             linetype = "dotted") +
  geom_vline(aes(xintercept = 0),
             color = "red",
             linetype = "dotted") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0, 0.5, 1))+
  scale_color_manual(values = c("grey60", "#7AD151FF", "#2A788EFF", "#414487FF")) +
  scale_fill_manual(values = c("grey60", "#7AD151FF",  "#2A788EFF", "#414487FF")) +
  ylab("density") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "top",
        axis.title.x = element_blank(),
        legend.text = element_text(size = 16)) 

}  
  
