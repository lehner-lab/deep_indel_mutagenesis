calculate_meansub_pos <- function(x){
  
  meansub_pos<-data.frame(aggregate(scaled_variants_aPCA[scaled_variants_aPCA$mut_type=="substitutions" &
                                                                     scaled_variants_aPCA$domain == x,]$scaled_fitness,
                                              by=list(scaled_variants_aPCA[scaled_variants_aPCA$mut_type=="substitutions" &
                                                                             scaled_variants_aPCA$domain == x,]$Pos), FUN = mean))
  names(meansub_pos)<-c("Pos", "scaled_fitness_mean")
  
  ## add range of effects by calculating sd per position
  sdsub_pos<-data.frame(aggregate(scaled_variants_aPCA[scaled_variants_aPCA$mut_type=="substitutions" &
                                                                   scaled_variants_aPCA$domain == x,]$scaled_fitness,
                                            by=list(scaled_variants_aPCA[scaled_variants_aPCA$mut_type=="substitutions" &
                                                                           scaled_variants_aPCA$domain == x,]$Pos), FUN = sd))
  names(sdsub_pos)<-c("Pos", "scaled_fitness_sd")
  
  # merge to create one df
  meansub_pos<-merge(meansub_pos,
                     sdsub_pos, by="Pos")
  
  return(meansub_pos = meansub_pos)
  
}
