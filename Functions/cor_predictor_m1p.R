cor_predictor_m1p <- function(){
  
  df<-c()
  for (i in 1:length(file_list)){
    if (nrow(file_list[[i]])>3){
      temp<-cor.test(file_list[[i]]$ddG,
                     file_list[[i]]$model1p,
                     method="pearson")
      
      temp2<-cor.test(file_list[[i]]$ddG,
                      file_list[[i]]$model1p,
                      method="spearman")
      
      temp3<-data.frame(domain=unique(file_list[[i]]$domain),
                        Pearson=temp$estimate,
                        Pearson_pvalue=temp$p.value,
                        Spearman=temp2$estimate,
                        Spearman_pvalue=temp2$p.value)
      df <- df %>% rbind(temp3)
    }
  }
  
  return(list(df))
}
