scatter_figure_3_rSASA_subs <- function(x){
  
  ggplot(x,
         aes(x=rSASA,
             y=scaled_fitness_mean,
             color=domain))+
    geom_point()+
    theme_classic()+
    ylab("abundance (mean sub/residue)")+
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
