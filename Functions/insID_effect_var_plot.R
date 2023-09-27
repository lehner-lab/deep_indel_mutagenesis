insID_effect_var_plot <- function(){
  
  ggplot(data = summary_df, aes(x = Pos, y = (min_fitness + max_fitness) / 2)) +
    geom_linerange(aes(ymin = min_fitness, ymax = max_fitness),
                   position = position_dodge(width = 0.2)) +
    geom_point(data = df[df$insID != "A", ],
               aes(y = scaled_fitness),
               color = "grey") +
    geom_point(data = df[df$insID == "A", ],
               aes(y = scaled_fitness),
               color = "#2A788EFF",
               position = position_dodge(width = 0.2)) +
    geom_text(aes(x=label_df$Pos,
                  y=label_df$max_fitness+0.1,
                  label=label_df$max_insID))+
    geom_text(aes(x=label_df_min$Pos,
                  y=label_df_min$min_fitness-0.1,
                  label=label_df_min$min_insID))+
    labs(x = "position in domain", y = "protein abundance") +
    scale_x_continuous(breaks = summary_df$Pos, labels = summary_df$Pos) +
    theme_minimal()+
    theme(axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=16),
          axis.title.x  = element_text(size=16),
          axis.title.y = element_blank())
  
}
