
plot_validation <- function(domain){
  
    df = cor_tsuboyamaVStopolska_realigned[cor_tsuboyamaVStopolska_realigned$domain == domain,]
    
    ggplot()+
      geom_point(data=df,
             aes(x=deltaG_ins,
                 y=scaled_fitness_ins,
                 color="ins"))+
      geom_point(data=df,
             aes(x=deltaG_dels,
                 y=scaled_fitness_dels,
                 color="del"))+
      geom_text_repel(data=df,
                  aes(x=deltaG_ins,
                      y=scaled_fitness_dels,
                      color="ins",
                      label=Pos), size=5)+
      geom_text_repel(data=df,
                  aes(x=deltaG_dels,
                      y=scaled_fitness_dels,
                      color="ins",
                      label=Pos), size=5)+
      scale_color_manual(values=c("#5DC863FF", "#21908CFF", "#3B528BFF", "#f9a242ff"))+
      theme_classic()+
      ylab("Topolska et al., protein abundance")+
      xlab("Tsuboyama et al., ddG")+
      theme(axis.title = element_blank(),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 18, margin = margin(b = 2)),
        legend.position = "none")
}
