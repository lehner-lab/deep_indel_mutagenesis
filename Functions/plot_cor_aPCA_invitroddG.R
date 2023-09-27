plot_cor_aPCA_invitroddG<-function(x,y){
  ggplot(data=x,
         aes(x=ddG,
             y=growthrate,
             color=domain))+
    geom_pointrange(data=x,
                    aes(x=ddG,
                        y=growthrate,
                        color=domain,
                        ymin=growthrate-growthrate_sigma,
                        ymax=growthrate+growthrate_sigma))+
    geom_point()+
    #geom_text_repel(aes(size=2))+
    theme_classic()+
    xlab("mutant ddG (in vitro)")+
    ylab("abundance (growthrate)")+
    scale_color_manual(values=y)+
    theme(axis.title.y = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.text.y = element_text(size = 20, 
                                     vjust = 0.5, 
                                     hjust = 0.5, 
                                     angle = 0),
          axis.text.x = element_text(size = 20, 
                                     vjust = 0.5, 
                                     hjust = 0.5, 
                                     angle = 0),
          legend.position = "none",
          legend.text = element_text(size=26),
          legend.title = element_blank())
}
