load_coefficents_data_dels <- function(directory_path){
  
  directory_path <- directory_path
  
  # Get a list of all .rds files in the directory
  file_list <- list.files(directory_path, pattern = "\\.rds$", full.names = TRUE)
  
  
  list_coef<-c()
  for (i in 1:length(file_list)){
    pdb_name<-substr(file_list[[i]], 113, 116) ## adjust this number to the pdb name in the file_list path
    tmp<-readRDS(file_list[[i]])
    
    
    r_best_coef_m2_noin=tmp[[1]]
    r_best_coef_m3_noin=tmp[[2]]
    r_best_coef_m3_inter=tmp[[3]]
    r_best_coef_m4_noin=tmp[[4]]
    r_best_coef_m4_inter=tmp[[5]]
    r_best_coef_m5_noin=tmp[[6]]
    r_best_coef_m5_inter=tmp[[7]]
    r_best_coef_m5p_noin=tmp[[8]]
    r_best_coef_m5p_inter=tmp[[9]]
    r_best_coef_m6_noin=tmp[[10]]
    r_best_coef_m7_noin=tmp[[11]]
    r_best_coef_m7_inter=tmp[[12]]
    r_best_coef_m7p_noin=tmp[[13]]
    r_best_coef_m7p_inter=tmp[[14]]
    

    tmp2<-list(r_best_coef_m2_noin,
               r_best_coef_m3_noin,
               r_best_coef_m3_inter,
               r_best_coef_m4_noin,
               r_best_coef_m4_inter,
               r_best_coef_m5_noin,
               r_best_coef_m5_inter,
               r_best_coef_m5p_noin,
               r_best_coef_m5p_inter,
               r_best_coef_m6_noin,
               r_best_coef_m7_noin,
               r_best_coef_m7_inter,
               r_best_coef_m7p_noin,
               r_best_coef_m7p_inter)
    
    for (j in 1:length(tmp2)){
      coef_value<-as.data.frame(as.matrix(tmp2[[j]]))
      coef_name<-rownames(tmp2[[j]])
      
      tmp3<-data.frame(coef_name=coef_name,
                       coef_value=coef_value)
      rownames(tmp3)<-1:nrow(tmp3)
      
      tmp2[[j]]<-tmp3
    }
    
    models<-c("m2", "m3", "m3_inter", "m4", "m4_inter", "m5", "m5_inter", "m5p", "m5p_inter", "m6", "m7", "m7_inter", "m7p", "m7p_inter")
    
    for (k in 1:length(tmp2)){
      tmp2[[k]]$model<-""
      tmp2[[k]]$model<-models[[k]]
    }
    
    tmp<-bind_rows(tmp2)
    tmp$pdb_name<-pdb_name
    
    list_coef <- list_coef %>% bind_rows(tmp)
    
  }
  
  list_coef<-split(list_coef, f=list_coef$model)
  
  return(list_coef =list_coef)
  
}


