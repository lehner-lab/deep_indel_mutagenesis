plot_density_single_indels<-function(x){
  ggplot(x, 
         aes(x = scaled_fitness, y = mut_type))+
    geom_density_ridges(aes(color=mut_type,
                            scale = 1.6,
                            fill=mut_type,
                            alpha="#FF0000A0"))+
    geom_vline(xintercept=1,
               linetype="dotted")+
    geom_vline(xintercept=0, 
               color="red", 
               linetype="dotted")+
    scale_y_discrete(expand = c(0, 0))+
    scale_x_continuous(breaks=seq(-1,2,0.5))+
    scale_color_manual(values=c("deletions" = "#7AD151FF", "insertions" = "#2A788EFF", "substitutions" = "#414487FF"))+
    scale_fill_manual(values=c("deletions" = "#7AD151FF", "insertions" = "#2A788EFF", "substitutions" = "#414487FF"))+
    ylab("density")+
    xlab("protein abundance")+
    theme_classic()+
    theme(axis.title.y = element_blank(),
          axis.text.y =  element_blank(),
          axis.text.x = element_text(size=20),
          legend.title = element_blank(),
          legend.position = "none",
          axis.title.x = element_text(size=22),
          legend.text = element_text(size=20))+
    guides(alpha=FALSE,
           colour = guide_legend(reverse=T),
           fill = guide_legend(reverse=T))
  
  
}

