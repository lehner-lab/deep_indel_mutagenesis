scatter_for_fig2<-function(d,x,y,s){

ggplot(data=d,
       aes(x=x,
           y=y,
           label=Pos))+
  geom_errorbar(data=d,
                aes(x=x,
                    y=y,
                    ymin=y-s,
                    ymax=y+s, color="grey60"))+
  geom_text(aes(label=Pos,
                size=60,
                fontface="bold"))+
  theme_classic()+
  scale_color_manual(values=c("gray60","black"))+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 22, 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0),
        legend.position = "none")

}
