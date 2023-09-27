scatter_for_fig2_2 <- function(df1, df2){
  
  variation_data_mut_subs <- df1 %>%
    group_by(Mut) %>%
    summarize(
      mean_sub_pos = mean(scaled_fitness)
    )

    colnames(df2)[4] <- "Mut"

  variation_data_mut_ins <- df2 %>%
    group_by(Mut) %>%
    summarize(
      mean_ins_pos = mean(scaled_fitness)
    )
  
  merged_df <- merge(variation_data_mut_subs,
                     variation_data_mut_ins,
                     by="Mut")
  
  ggplot(data=merged_df,
       aes(x=mean_sub_pos,
           y=mean_ins_pos))+
  geom_text(aes(label=Mut,
                size=70,
                fontface="bold"))+
  theme_classic()+
  ylab("mean ins/mutation")+
  xlab("mean sub/mutation")+
  scale_color_manual(values=c("gray60","black"))+
  theme(axis.title.y = element_text(size=22),
        axis.title.x = element_text(size=22),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18, 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0),
        legend.position = "none")
}
