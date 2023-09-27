scatter_figure_3_rSASA_indels <- function(x){
  
  ggplot(x,
         aes(x=rSASA,
             y=scaled_fitness,
             color=domain))+
    geom_pointrange(data=x,
                    aes(x=rSASA,
                        y=scaled_fitness,
                        ymax=scaled_fitness+scaled_sigma,
                        ymin=scaled_fitness-scaled_sigma,
                        color=domain))+
    theme_classic()+
    xlab("rSASA")+
    theme(axis.title.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.text.x = element_text(size = 18, 
                                     vjust = 0.5, 
                                     hjust = 0.5, 
                                     angle = 0),
          legend.title = element_blank(),
          legend.position ="top",
          legend.text = element_text(size=18))+
    scale_color_manual(values=c("grey20", "grey80"))
  
}
