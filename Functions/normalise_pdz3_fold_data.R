normalise_pdz3_fold_data <- function() {
  
  ##extract position of mutation, add Pos  and length column
  pdz3_fold_synonymous$diff<-as.character(pdz3_fold_synonymous$diff)
  pdz3_fold_synonymous$diff<-str_first_number(pdz3_fold_synonymous$diff)
  pdz3_fold_synonymous$diff<-as.numeric(pdz3_fold_synonymous$diff)
  pdz3_fold_synonymous$Pos<-as.numeric(pdz3_fold_synonymous$diff)
  pdz3_fold_synonymous$lenght<-str_length(pdz3_fold_synonymous$aa_seq)
  
  
  ##extract positions as numeric into the Pos column for all df
  pdz3_fold_singleDEL$Pos<-as.numeric(pdz3_fold_singleDEL$diff)
  
  pdz3_fold_doubleDEL$diff<-as.character(pdz3_fold_doubleDEL$diff)
  pdz3_fold_doubleDEL$Pos<-str_first_number(pdz3_fold_doubleDEL$diff)
  pdz3_fold_doubleDEL$Pos<-as.numeric(pdz3_fold_doubleDEL$Pos)
  
  pdz3_fold_tripleDEL$diff<-as.character(pdz3_fold_tripleDEL$diff)
  pdz3_fold_tripleDEL$Pos<-str_first_number(pdz3_fold_tripleDEL$diff)
  pdz3_fold_tripleDEL$Pos<-as.numeric(pdz3_fold_tripleDEL$Pos)
  
  pdz3_fold_allins$Pos<-as.numeric(pdz3_fold_allins$diff)
  pdz3_fold_CX$Pos<-as.numeric(pdz3_fold_CX$diff)
  
  pdz3_fold_CXX$diff<-as.character(pdz3_fold_CXX$diff)
  pdz3_fold_CXX$Pos<-str_first_number(pdz3_fold_CXX$diff)
  pdz3_fold_CXX$Pos<-as.numeric(pdz3_fold_CXX$Pos)
  
  pdz3_fold_CXXX$diff<-as.character(pdz3_fold_CXXX$diff)
  pdz3_fold_CXXX$Pos<-str_first_number(pdz3_fold_CXXX$diff)
  pdz3_fold_CXXX$Pos<-as.numeric(pdz3_fold_CXXX$Pos)
  
  pdz3_fold_Delsub$diff<-as.character(pdz3_fold_Delsub$diff)
  pdz3_fold_Delsub$Pos<-str_first_number(pdz3_fold_Delsub$diff)
  
  ### find the aa position of the synonymous mutation
  pdz3_fold_synonymous$Pos <- as.character(floor(pdz3_fold_synonymous$diff / 3))
  
  
  ## add mutation type
  pdz3_fold_singleDEL$type<-"singleDEL"
  pdz3_fold_doubleDEL$type<-"doubleDEL"
  pdz3_fold_tripleDEL$type<-"tripleDEL"
  pdz3_fold_synonymous$type<-"synonymous"
  pdz3_fold_CX$type<-"singleINS"
  pdz3_fold_CXX$type<-"doubleINS"
  pdz3_fold_CXXX$type<-"tripleINS"
  pdz3_fold_Delsub$type<-"delSub"    
  pdz3_fold_allins$type<-"allins"
  
  ## data normalisation
  # normalise the data to mode of dead peak of the bimodal distribution of variants
  
  #bind all data
  test<-rbind(pdz3_fold_synonymous,
              pdz3_fold_allins,
              pdz3_fold_singles,
              pdz3_fold_CX,
              pdz3_fold_CXX,
              pdz3_fold_CNN,
              pdz3_fold_insAA,
              pdz3_fold_insCC,
              pdz3_fold_CXXX,
              pdz3_fold_CNNN,
              pdz3_fold_singleDEL,
              pdz3_fold_doubleDEL, 
              pdz3_fold_tripleDEL,
              fill=TRUE)
  
  ##find mode of the lower peak in the bimodal distribution --> referred to as STOPs
  STOPs<-min(mlv(na.omit(test[test$fitness<0,]$fitness), method="naive"))
  
  ## substract the mode from all variants. 
  pdz3_fold_CX$norm_fitness<-pdz3_fold_CX$fitness-STOPs
  pdz3_fold_CXX$norm_fitness<-pdz3_fold_CXX$fitness-STOPs
  pdz3_fold_CXXX$norm_fitness<-pdz3_fold_CXXX$fitness-STOPs
  pdz3_fold_singleDEL$norm_fitness<-pdz3_fold_singleDEL$fitness-STOPs
  pdz3_fold_doubleDEL$norm_fitness<-pdz3_fold_doubleDEL$fitness-STOPs
  pdz3_fold_tripleDEL$norm_fitness<-pdz3_fold_tripleDEL$fitness-STOPs
  pdz3_fold_Delsub$norm_fitness<-pdz3_fold_Delsub$fitness-STOPs
  pdz3_fold_synonymous$norm_fitness<-pdz3_fold_synonymous$fitness-STOPs
  pdz3_fold_insAA$norm_fitness<-pdz3_fold_insAA$fitness-STOPs
  pdz3_fold_insCC$norm_fitness<-pdz3_fold_insCC$fitness-STOPs
  pdz3_fold_CNN$norm_fitness<-pdz3_fold_CNN$fitness-STOPs
  pdz3_fold_CNNN$norm_fitness<-pdz3_fold_CNNN$fitness-STOPs
  pdz3_fold_allins$norm_fitness<-pdz3_fold_allins$fitness-STOPs
  pdz3_fold_singles$norm_fitness<-pdz3_fold_singles$fitness-STOPs
  
  ## find weighted mean of the synonymous variants
  syns<-sum(pdz3_fold_synonymous$norm_fitness/(pdz3_fold_synonymous$sigma^2))/sum(1/(pdz3_fold_synonymous$sigma^2))
  
  ##normalise by dividing variants with the weighted mean of the synonymous variants
  pdz3_fold_synonymous$scaled_fitness<-pdz3_fold_synonymous$norm_fitness/syns
  pdz3_fold_CX$scaled_fitness<-pdz3_fold_CX$norm_fitness/syns
  pdz3_fold_CXX$scaled_fitness<-pdz3_fold_CXX$norm_fitness/syns
  pdz3_fold_CXXX$scaled_fitness<-pdz3_fold_CXXX$norm_fitness/syns
  pdz3_fold_singleDEL$scaled_fitness<-pdz3_fold_singleDEL$norm_fitness/syns
  pdz3_fold_doubleDEL$scaled_fitness<-pdz3_fold_doubleDEL$norm_fitness/syns
  pdz3_fold_tripleDEL$scaled_fitness<-pdz3_fold_tripleDEL$norm_fitness/syns
  pdz3_fold_Delsub$scaled_fitness<-pdz3_fold_Delsub$norm_fitness/syns
  pdz3_fold_insAA$scaled_fitness<-pdz3_fold_insAA$norm_fitness/syns
  pdz3_fold_insCC$scaled_fitness<-pdz3_fold_insCC$norm_fitness/syns
  pdz3_fold_CNN$scaled_fitness<-pdz3_fold_CNN$norm_fitness/syns
  pdz3_fold_CNNN$scaled_fitness<-pdz3_fold_CNNN$norm_fitness/syns
  pdz3_fold_allins$scaled_fitness<-pdz3_fold_allins$norm_fitness/syns
  pdz3_fold_singles$scaled_fitness<-pdz3_fold_singles$norm_fitness/syns
  
  ## scale the error by dividing by the same factor: dividing by the weighted mean of the synonmous fitness
  pdz3_fold_synonymous$scaled_sigma<-pdz3_fold_synonymous$sigma/(syns)
  pdz3_fold_CX$scaled_sigma<-pdz3_fold_CX$sigma/(syns)
  pdz3_fold_CXX$scaled_sigma<-pdz3_fold_CXX$sigma/(syns)
  pdz3_fold_CXXX$scaled_sigma<-pdz3_fold_CXXX$sigma/(syns)
  pdz3_fold_singleDEL$scaled_sigma<-pdz3_fold_singleDEL$sigma/(syns)
  pdz3_fold_doubleDEL$scaled_sigma<-pdz3_fold_doubleDEL$sigma/(syns)
  pdz3_fold_tripleDEL$scaled_sigma<-pdz3_fold_tripleDEL$sigma/(syns)
  pdz3_fold_Delsub$scaled_sigma<-pdz3_fold_Delsub$sigma/(syns)
  pdz3_fold_insAA$scaled_sigma<-pdz3_fold_insAA$sigma/(syns)
  pdz3_fold_insCC$scaled_sigma<-pdz3_fold_insCC$sigma/(syns)
  pdz3_fold_CNN$scaled_sigma<-pdz3_fold_CNN$sigma/(syns)
  pdz3_fold_CNNN$scaled_sigma<-pdz3_fold_CNNN$sigma/(syns)
  pdz3_fold_allins$scaled_sigma<-pdz3_fold_allins$sigma/(syns)
  pdz3_fold_singles$scaled_sigma<-pdz3_fold_singles$sigma/(syns)
  

  ##fix mutation positions that are called wrong manually. 
  # this is just an artifact of the function used to call the positions of difference and differing aa sequence
  pdz3_fold_doubleDEL[pdz3_fold_doubleDEL$diff == "c(2, 4)",]$Pos<-3
  pdz3_fold_doubleDEL[pdz3_fold_doubleDEL$diff == "c(19, 21)",]$Pos<-20
  pdz3_fold_doubleDEL[pdz3_fold_doubleDEL$diff == "c(34, 36)",]$Pos<-35
  
  pdz3_fold_tripleDEL[pdz3_fold_tripleDEL$diff == "c(2, 4, 5)",]$Pos<-3
  pdz3_fold_tripleDEL[pdz3_fold_tripleDEL$diff == "c(4, 5, 7)",]$Pos<-5
  pdz3_fold_tripleDEL[pdz3_fold_tripleDEL$diff == "c(12, 13, 15)",]$Pos<-13
  pdz3_fold_tripleDEL[pdz3_fold_tripleDEL$diff == "c(19, 21, 22)",]$Pos<-20
  pdz3_fold_tripleDEL[pdz3_fold_tripleDEL$diff == "c(23, 24, 26)",]$Pos<-24
  pdz3_fold_tripleDEL[pdz3_fold_tripleDEL$diff == "c(26, 27, 29)",]$Pos<-27
  pdz3_fold_tripleDEL[pdz3_fold_tripleDEL$diff == "c(34, 36, 37)",]$Pos<-35
  
  ## add missing variants to the df
  # these are variants that have the same aa_seq although the mutation is technically in different positons. 
  #singleDEL:2=3,19=20,34=35
  pdz3_fold_singleDEL_missing<-pdz3_fold_singleDEL[c(2,8,14),]
  pdz3_fold_singleDEL_missing$Pos<-c(3,20,35)
  pdz3_fold_singleDEL<-rbind(pdz3_fold_singleDEL,
                             pdz3_fold_singleDEL_missing)
  
  #doubleDEL: 4=5, 12=13, 23=24, 26=27
  pdz3_fold_doubleDEL_missing<-pdz3_fold_doubleDEL[c(1,40,10,11),]
  pdz3_fold_doubleDEL_missing$Pos<-c(5,13,24,27)
  pdz3_fold_doubleDEL<-rbind(pdz3_fold_doubleDEL,
                             pdz3_fold_doubleDEL_missing)
  
  #tripleDEL:9=10, 20=21, 20=22, 27=28, 27=29
  pdz3_fold_tripleDEL_missing<-pdz3_fold_tripleDEL[c(5,11,11,14,14),]
  pdz3_fold_tripleDEL_missing$Pos<-c(10,21,22,28,29)
  pdz3_fold_tripleDEL<-rbind(pdz3_fold_tripleDEL,
                             pdz3_fold_tripleDEL_missing)
  #CX 2=3,19=20,34=35
  pdz3_fold_CX_missing<-pdz3_fold_CX[c(2,43,37),]
  pdz3_fold_CX_missing$Pos<-c(3,20,35)
  pdz3_fold_CX<-rbind(pdz3_fold_CX,
                      pdz3_fold_CX_missing)
  
  

  #### test if the change in fitness is signifcant
  ## we need to test for significance for values > and < than 1.  
  ## We calculate a z-stat and then do a 2-tailed test
  # we do mulitple testing correction using bonferroni
  
  
  # CX
  zstats=(pdz3_fold_CX$scaled_fitness-1)/pdz3_fold_CX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  pdz3_fold_CX$zstat=zstats
  pdz3_fold_CX$pvalue=pvals_man
  pdz3_fold_CX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  pdz3_fold_CX$fdrsig=pdz3_fold_CX$FDR<=.05
  pdz3_fold_CX<-pdz3_fold_CX %>% mutate(significant=
                                          ifelse(fdrsig == "TRUE", "*",NA))
  
  
  ##CXX
  zstats=(pdz3_fold_CXX$scaled_fitness-1)/pdz3_fold_CXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  pdz3_fold_CXX$zstat=zstats
  pdz3_fold_CXX$pvalue=pvals_man
  pdz3_fold_CXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  pdz3_fold_CXX$fdrsig=pdz3_fold_CXX$FDR<=.05
  pdz3_fold_CXX<-pdz3_fold_CXX %>% mutate(significant=
                                            ifelse(fdrsig == "TRUE", "*",NA))
  ##CXXX
  zstats=(pdz3_fold_CXXX$scaled_fitness-1)/pdz3_fold_CXXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  pdz3_fold_CXXX$zstat=zstats
  pdz3_fold_CXXX$pvalue=pvals_man
  pdz3_fold_CXXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  pdz3_fold_CXXX$fdrsig=pdz3_fold_CXXX$FDR<=.05
  pdz3_fold_CXXX<-pdz3_fold_CXXX %>% mutate(significant=
                                              ifelse(fdrsig == "TRUE", "*",NA))
  #singleDEL
  zstats=(pdz3_fold_singleDEL$scaled_fitness-1)/pdz3_fold_singleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  pdz3_fold_singleDEL$zstat=zstats
  pdz3_fold_singleDEL$pvalue=pvals_man
  pdz3_fold_singleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  pdz3_fold_singleDEL$fdrsig=pdz3_fold_singleDEL$FDR<=.05
  pdz3_fold_singleDEL<-pdz3_fold_singleDEL %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #doubleDEL
  zstats=(pdz3_fold_doubleDEL$scaled_fitness-1)/pdz3_fold_doubleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  pdz3_fold_doubleDEL$zstat=zstats
  pdz3_fold_doubleDEL$pvalue=pvals_man
  pdz3_fold_doubleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  pdz3_fold_doubleDEL$fdrsig=pdz3_fold_doubleDEL$FDR<=.05
  pdz3_fold_doubleDEL<-pdz3_fold_doubleDEL %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #tripleDEL
  zstats=(pdz3_fold_tripleDEL$scaled_fitness-1)/pdz3_fold_tripleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  pdz3_fold_tripleDEL$zstat=zstats
  pdz3_fold_tripleDEL$pvalue=pvals_man
  pdz3_fold_tripleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  pdz3_fold_tripleDEL$fdrsig=pdz3_fold_tripleDEL$FDR<=.05
  pdz3_fold_tripleDEL<-pdz3_fold_tripleDEL %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #DelSub
  zstats=(pdz3_fold_Delsub$scaled_fitness-1)/pdz3_fold_Delsub$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  pdz3_fold_Delsub$zstat=zstats
  pdz3_fold_Delsub$pvalue=pvals_man
  pdz3_fold_Delsub[,FDR:=p.adjust(pvalue,method="bonferroni")]
  pdz3_fold_Delsub$fdrsig=pdz3_fold_Delsub$FDR<=.05
  pdz3_fold_Delsub<-pdz3_fold_Delsub %>% mutate(significant=
                                                  ifelse(fdrsig == "TRUE", "*",NA))
  
  pdz3_fold_Delsub<-pdz3_fold_Delsub[pdz3_fold_Delsub$STOP==FALSE,]
  
  
  ##some of the DelSubs can cause 2 different substitutions at the same positions. We need to seperate them out into 2 dataframes. 
  pdz3_fold_Delsub_2<-pdz3_fold_Delsub[duplicated(pdz3_fold_Delsub$Pos),]
  pdz3_fold_Delsub<-pdz3_fold_Delsub[!duplicated(pdz3_fold_Delsub$Pos),]
  pdz3_fold_Delsub_2$type<-"delSub_1"
  pdz3_fold_Delsub$type<-"delSub_2"
  
  #### CNN
  zstats=(pdz3_fold_CNN$scaled_fitness-1)/pdz3_fold_CNN$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  pdz3_fold_CNN$zstat=zstats
  pdz3_fold_CNN$pvalue=pvals_man
  pdz3_fold_CNN[,FDR:=p.adjust(pvalue,method="bonferroni")]
  pdz3_fold_CNN$fdrsig=pdz3_fold_CNN$FDR<=.05
  pdz3_fold_CNN<-pdz3_fold_CNN %>% mutate(significant=
                                            ifelse(fdrsig == "TRUE", "*",NA))
  
  pdz3_fold_CNN$Pos<-str_first_number(pdz3_fold_CNN$diff)
  
  #### insAA
  zstats=(pdz3_fold_insAA$scaled_fitness-1)/pdz3_fold_insAA$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  pdz3_fold_insAA$zstat=zstats
  pdz3_fold_insAA$pvalue=pvals_man
  pdz3_fold_insAA[,FDR:=p.adjust(pvalue,method="bonferroni")]
  pdz3_fold_insAA$fdrsig=pdz3_fold_insAA$FDR<=.05
  pdz3_fold_insAA<-pdz3_fold_insAA %>% mutate(significant=
                                                ifelse(fdrsig == "TRUE", "*",NA))
  
  pdz3_fold_insAA$Pos<-str_first_number(pdz3_fold_insAA$diff)
  pdz3_fold_insAA$type<-"doubleAA"
  
  ### insCC
  zstats=(pdz3_fold_insCC$scaled_fitness-1)/pdz3_fold_insCC$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  pdz3_fold_insCC$zstat=zstats
  pdz3_fold_insCC$pvalue=pvals_man
  pdz3_fold_insCC[,FDR:=p.adjust(pvalue,method="bonferroni")]
  pdz3_fold_insCC$fdrsig=pdz3_fold_insCC$FDR<=.05
  pdz3_fold_insCC<-pdz3_fold_insCC %>% mutate(significant=
                                                ifelse(fdrsig == "TRUE", "*",NA))
  
  pdz3_fold_insCC$Pos<-str_first_number(pdz3_fold_insCC$diff)
  pdz3_fold_insCC$type<-"doubleCC"
  
  ### CNNN
  zstats=(pdz3_fold_CNNN$scaled_fitness-1)/pdz3_fold_CNNN$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  pdz3_fold_CNNN$zstat=zstats
  pdz3_fold_CNNN$pvalue=pvals_man
  pdz3_fold_CNNN[,FDR:=p.adjust(pvalue,method="bonferroni")]
  pdz3_fold_CNNN$fdrsig=pdz3_fold_CNNN$FDR<=.05
  pdz3_fold_CNNN<-pdz3_fold_CNNN %>% mutate(significant=
                                              ifelse(fdrsig == "TRUE", "*",NA))
  
  pdz3_fold_CNNN$Pos<-str_first_number(pdz3_fold_CNNN$diff)
  
  ### allins
  zstats=(pdz3_fold_allins$scaled_fitness-1)/pdz3_fold_allins$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  pdz3_fold_allins$zstat=zstats
  pdz3_fold_allins$pvalue=pvals_man
  pdz3_fold_allins[,FDR:=p.adjust(pvalue,method="bonferroni")]
  pdz3_fold_allins$fdrsig=pdz3_fold_allins$FDR<=.05
  pdz3_fold_allins<-pdz3_fold_allins %>% mutate(significant=
                                                  ifelse(fdrsig == "TRUE", "*",NA))
  ### singles
  zstats=(pdz3_fold_singles$scaled_fitness-1)/pdz3_fold_singles$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  pdz3_fold_singles$zstat=zstats
  pdz3_fold_singles$pvalue=pvals_man
  pdz3_fold_singles[,FDR:=p.adjust(pvalue,method="bonferroni")]
  pdz3_fold_singles$fdrsig=pdz3_fold_singles$FDR<=.05
  pdz3_fold_singles<-pdz3_fold_singles %>% mutate(significant=
                                                    ifelse(fdrsig == "TRUE", "*",NA))
  pdz3_fold_singles$type<-"substitution"
  
  
  
  ## bind all variants. 
  pdz3_fold_allvariants<-rbind(pdz3_fold_synonymous[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
                               pdz3_fold_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               pdz3_fold_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               pdz3_fold_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               pdz3_fold_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               pdz3_fold_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")], 
                               pdz3_fold_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               pdz3_fold_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               pdz3_fold_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               pdz3_fold_singles[pdz3_fold_singles$STOP=="FALSE",c("aa_seq", "Pos", "Mut", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               fill=TRUE)
  
  
  return(list(pdz3_fold_allvariants = pdz3_fold_allvariants, 
              pdz3_fold_CX = pdz3_fold_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              pdz3_fold_CXX = pdz3_fold_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              pdz3_fold_CXXX = pdz3_fold_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              pdz3_fold_CNN = pdz3_fold_CNN[,c("aa_seq", "Pos","insID", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              pdz3_fold_CNNN = pdz3_fold_CNNN[,c("aa_seq", "Pos", "insID","lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              pdz3_fold_insCC = pdz3_fold_insCC[,c("aa_seq", "Pos","insID", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              pdz3_fold_insAA = pdz3_fold_insAA[,c("aa_seq", "Pos","insID", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              pdz3_fold_allins = pdz3_fold_allins[,c("aa_seq", "Pos","insID", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              pdz3_fold_synonymous = pdz3_fold_synonymous[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
              pdz3_fold_singles = pdz3_fold_singles[,c("aa_seq", "Pos", "Mut", "type", "scaled_fitness", "scaled_sigma", "significant")],
              pdz3_fold_Delsub = pdz3_fold_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              pdz3_fold_Delsub_2 = pdz3_fold_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              pdz3_fold_singleDEL = pdz3_fold_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              pdz3_fold_doubleDEL = pdz3_fold_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              pdz3_fold_tripleDEL = pdz3_fold_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")]))
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
