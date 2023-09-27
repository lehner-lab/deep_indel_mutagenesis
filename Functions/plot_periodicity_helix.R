plot_periodity_helix <- function(ddG, pos, y, counts){
  
  ## find mean
  ddG_mean<-aggregate(ddG, by=list(pos),FUN=mean)
  ddG_mean<-as.data.frame(ddG_mean)
  colnames(ddG_mean)[1]<-"align_to_center"
  colnames(ddG_mean)[2]<-"mean_ddG"

ggplot(data = df, 
       aes(x = reorder(as.character(align_to_center), 
                       as.numeric(align_to_center)), 
           y = ddG)) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.02) +
  scale_y_reverse() +
  geom_point(data = ddG_mean,
             aes(x=reorder(as.character(align_to_center), 
                           as.numeric(align_to_center)),
                 y=mean_ddG,
                 color="red3"))+
  geom_line(data = ddG_mean,
            aes(x=reorder(as.character(align_to_center), 
                          as.numeric(align_to_center)),
                y=mean_ddG,
                color="red3",
                group=1))+
  geom_text(data = counts, 
            aes(x = reorder(as.character(align_to_center), as.numeric(align_to_center)), 
                y = -2, label = element_no_simple), 
            vjust = 1, size = 3)+
  theme_classic()+
  scale_color_manual(values=c("red3","red3"))+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=13),
        legend.position = "none")
  
}
