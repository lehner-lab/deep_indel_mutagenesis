make_lineplots_single_effects <- function(x,y,z){
  
    ggplot()+
      geom_rect(data=data_breaks,
                aes(xmin = start,
                    xmax = end,
                    ymin = - Inf,
                    ymax = Inf,
                    fill = colors),
                alpha = 0.5)+
      geom_point(data=x, 
                 aes(x=Pos,
                     y=scaled_fitness,
                     color="insertions"))+
      geom_line(data=x, 
                aes(x=Pos,
                    y=scaled_fitness,
                    color="insertions"))+
      geom_linerange(data=x, 
                     aes(x=Pos,
                         y=scaled_fitness,
                         ymin=scaled_fitness-scaled_sigma,
                         ymax=scaled_fitness+scaled_sigma,
                         color="insertions"))+
      geom_point(data=y, 
                 aes(x=Pos,
                     y=scaled_fitness,
                     color="deletions"))+
      geom_line(data=y, 
                aes(x=Pos,
                    y=scaled_fitness,
                    color="deletions"))+
      geom_linerange(data=y, 
                     aes(x=Pos,
                         y=scaled_fitness,
                         ymin=scaled_fitness-scaled_sigma,
                         ymax=scaled_fitness+scaled_sigma,
                         color="deletions"))+
      scale_color_manual(values=c("#7AD151FF","#2A788EFF"))+
      scale_fill_manual(values=c("red3","blue3","orange"))+
      xlab("position in domain")+
      ylab("protein abundance")+
      theme_classic()+
      scale_y_continuous(limits=c(-0.63,1.41))+
    scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 30, 35,40,45,50,55,60,65,70), 
                       labels = c(1, 5, 10, 15, 20, 25, 30, 35,40,45,50,55,60,65,70),
                       expand = c(0, 0),
                       limits = c(1, max(z$Pos)))+
      geom_hline(yintercept = 1,
                 color="black",
                 linetype="dotted")+
      geom_hline(yintercept = 0,
                 color="red",
                 linetype="dotted")+
    theme(axis.text.y = element_text(size = 20, 
                                     vjust = 0.5, 
                                     hjust = 0.5, 
                                     angle = 0),
          axis.text.x = element_text(size = 20, 
                                     vjust = 0.5, 
                                     hjust = 0.5, 
                                     angle = 0),
          legend.text = element_text(size=12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
    
}
