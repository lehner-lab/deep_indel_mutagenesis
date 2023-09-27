add_domain_name_and_join_dfs <- function() {
  
  O75400_PF01846_allvariants$domain <- "FBP11-FF1"
  P0A9X9_PF00313_allvariants$domain <- "CSPA-CSD"
  P01053_PF00280_allvariants$domain <- "CI2A-PIN1"
  P02417_PF01281_allvariants$domain <- "BL17-NTL9"
  P02640_PF02209_allvariants$domain <- "VIL1-HP"
  P32081_PF00313_allvariants$domain <- "CSPB-CSD"
  P61024_PF01111_allvariants$domain <- "CKS1"
  grb2_fold_allvariants$domain <- "GRB2-SH3"
  pdz3_fold_allvariants$domain <- "PSD95-PDZ3"
  
  grb2_fold_allins$domain <- "GRB2-SH3"
  grb2_fold_insAA$domain <- "GRB2-SH3"
  grb2_fold_insCC$domain <- "GRB2-SH3"
  grb2_fold_CNN$domain <- "GRB2-SH3"
  grb2_fold_CNNN$domain <- "GRB2-SH3"
  
  pdz3_fold_allins$domain <- "PSD95-PDZ3"
  pdz3_fold_insAA$domain <- "PSD95-PDZ3"
  pdz3_fold_insCC$domain <- "PSD95-PDZ3"
  pdz3_fold_CNN$domain <- "PSD95-PDZ3"
  pdz3_fold_CNNN$domain <- "PSD95-PDZ3"

  scaled_variants_aPCA <- rbind(O75400_PF01846_allvariants,
                                P0A9X9_PF00313_allvariants,
                                P01053_PF00280_allvariants,
                                P02417_PF01281_allvariants,
                                P02640_PF02209_allvariants,
                                P32081_PF00313_allvariants,
                                P61024_PF01111_allvariants,
                                grb2_fold_allvariants,
                                pdz3_fold_allvariants,
                                grb2_fold_allins,
                                pdz3_fold_allins,
                                grb2_fold_insAA,
                                grb2_fold_insCC,
                                pdz3_fold_insAA,
                                pdz3_fold_insCC,
                                grb2_fold_CNN,
                                grb2_fold_CNNN,
                                pdz3_fold_CNN,
                                pdz3_fold_CNNN,
                                fill=T)
  
  scaled_variants_aPCA <- scaled_variants_aPCA %>% mutate(mut_type=
                                                        ifelse(type == "synonymous", "synonymous",
                                                               ifelse(type == "substitution", "substitutions",
                                                                      ifelse(type == "STOP", "STOP",
                                                                             ifelse(type == "singleINS"| type == "doubleINS" | type == "tripleINS", "insertions",
                                                                                    ifelse(type == "singleDEL"| type == "doubleDEL" | type == "tripleDEL", "deletions",
                                                                                           ifelse(type == "delSub_1"| type == "delSub_2"| type == "delSub", "delSub",NA)))))))
  
  # add wt
  for (i in 1:nrow(scaled_variants_aPCA)){
    if (scaled_variants_aPCA[i,]$Pos == 0){
      scaled_variants_aPCA[i,]$type <- "wt"
    }
  }
  
  
  return(list(scaled_variants_aPCA = scaled_variants_aPCA)) 
}
