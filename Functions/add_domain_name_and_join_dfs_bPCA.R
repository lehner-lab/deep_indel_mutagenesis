add_domain_name_and_join_dfs_bPCA <- function() {
  
 ## all GRB2-SH3
  grb2_bind_wildtype$type<-"wt"
  grb2_bPCA<-rbind(grb2_bind_CX,
                   grb2_bind_CXX,
                   grb2_bind_CNN,
                   grb2_bind_insAA,
                   grb2_bind_insCC,
                   grb2_bind_CXXX,
                   grb2_bind_CNNN,
                   grb2_bind_singles,
                   grb2_bind_allins,
                   grb2_bind_Delsub,
                   grb2_bind_Delsub_2,
                   grb2_bind_wildtype,
                   grb2_bind_singleDEL,
                   grb2_bind_doubleDEL,
                   grb2_bind_tripleDEL,
                   grb2_bind_synonymous,
                   fill=T)
  
  grb2_bPCA$domain <- "GRB2-SH3"
  
 ## all PSD95-PDZ3
  pdz3_bind_wildtype$type<-"wt"
  pdz3_bPCA<-rbind(pdz3_bind_CX,
                   pdz3_bind_CXX,
                   pdz3_bind_CNN,
                   pdz3_bind_insAA,
                   pdz3_bind_insCC,
                   pdz3_bind_CXXX,
                   pdz3_bind_CNNN,
                   pdz3_bind_singles,
                   pdz3_bind_allins,
                   pdz3_bind_Delsub,
                   pdz3_bind_Delsub_2,
                   pdz3_bind_wildtype,
                   pdz3_bind_singleDEL,
                   pdz3_bind_doubleDEL,
                   pdz3_bind_tripleDEL,
                   pdz3_bind_synonymous,
                   fill=T)
  
  pdz3_bPCA$domain <- "PSD95-PDZ3"
  
  ## bind together
  scaled_variants_bPCA <- rbind(grb2_bPCA,
                                pdz3_bPCA)
  
  scaled_variants_bPCA <- scaled_variants_bPCA %>% mutate(mut_type=
                                                            ifelse(type == "synonymous", "synonymous",
                                                                   ifelse(type == "wt", "wt",
                                                                          ifelse(type == "substitution", "substitutions",
                                                                                 ifelse(type == "STOP", "STOP",
                                                                                        ifelse(type == "singleINS"| type == "doubleINS" | type == "tripleINS", "insertions",
                                                                                               ifelse(type == "singleDEL"| type == "doubleDEL" | type == "tripleDEL", "deletions",
                                                                                                      ifelse(type == "delSub_1"| type == "delSub_2"| type == "delSub", "delSub",NA))))))))
  
  return(list(scaled_variants_bPCA = scaled_variants_bPCA[,c("Pos", "domain", "aa_seq", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant", "Mut", "insID", "mut_type")])) 
}
