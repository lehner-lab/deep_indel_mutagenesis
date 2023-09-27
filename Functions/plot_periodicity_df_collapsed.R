plot_periodity_df_collapsed <- function(df,min_interval,min,
                                        max_interval,max,factored_pos){
  
  ##### collapse the last residues
  # Collapse data for align_to_center in lower end
  df$align_to_center <- ifelse(df$align_to_center %in% min_interval,
    min,
    as.character(df$align_to_center)
  )
  
  # Collapse data for align_to_center in higer end
  df$align_to_center <- ifelse(df$align_to_center %in% max_interval,
                               max,
                               as.character(df$align_to_center)
  )
  
  
  # Create a new factor variable for align_to_center with desired order
  df$align_to_center <- factor(
    df$align_to_center,
    levels = factored_pos
  )
  
  return(list(df = df))
  
}
