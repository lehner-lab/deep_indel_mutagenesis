normalise_grb2_fold_data <- function() {
  
  ##extract position of mutation, add Pos  and length column
  grb2_fold_synonymous$diff<-as.character(grb2_fold_synonymous$diff)
  grb2_fold_synonymous$diff<-str_first_number(grb2_fold_synonymous$diff)
  grb2_fold_synonymous$diff<-as.numeric(grb2_fold_synonymous$diff)
  grb2_fold_synonymous$Pos<-as.numeric(grb2_fold_synonymous$diff)
  grb2_fold_synonymous$lenght<-str_length(grb2_fold_synonymous$aa_seq)
  
  
  ##extract positions as numeric into the Pos column for all df
  grb2_fold_singleDEL$Pos<-as.numeric(grb2_fold_singleDEL$diff)
  
  grb2_fold_doubleDEL$diff<-as.character(grb2_fold_doubleDEL$diff)
  grb2_fold_doubleDEL$Pos<-str_first_number(grb2_fold_doubleDEL$diff)
  grb2_fold_doubleDEL$Pos<-as.numeric(grb2_fold_doubleDEL$Pos)
  
  grb2_fold_tripleDEL$diff<-as.character(grb2_fold_tripleDEL$diff)
  grb2_fold_tripleDEL$Pos<-str_first_number(grb2_fold_tripleDEL$diff)
  grb2_fold_tripleDEL$Pos<-as.numeric(grb2_fold_tripleDEL$Pos)

  grb2_fold_allins$Pos<-as.numeric(grb2_fold_allins$diff)
  grb2_fold_CX$Pos<-as.numeric(grb2_fold_CX$diff)
  
  grb2_fold_CXX$diff<-as.character(grb2_fold_CXX$diff)
  grb2_fold_CXX$Pos<-str_first_number(grb2_fold_CXX$diff)
  grb2_fold_CXX$Pos<-as.numeric(grb2_fold_CXX$Pos)
  
  grb2_fold_CXXX$diff<-as.character(grb2_fold_CXXX$diff)
  grb2_fold_CXXX$Pos<-str_first_number(grb2_fold_CXXX$diff)
  grb2_fold_CXXX$Pos<-as.numeric(grb2_fold_CXXX$Pos)

  grb2_fold_Delsub$diff<-as.character(grb2_fold_Delsub$diff)
  grb2_fold_Delsub$Pos<-str_first_number(grb2_fold_Delsub$diff)
  
  ### find the aa position of the synonymous mutation
  grb2_fold_synonymous$Pos <- as.character(floor(grb2_fold_synonymous$diff / 3))
  
  ## add mutation type
  grb2_fold_singleDEL$type<-"singleDEL"
  grb2_fold_doubleDEL$type<-"doubleDEL"
  grb2_fold_tripleDEL$type<-"tripleDEL"
  grb2_fold_synonymous$type<-"synonymous"
  grb2_fold_CX$type<-"singleINS"
  grb2_fold_CXX$type<-"doubleINS"
  grb2_fold_CXXX$type<-"tripleINS"
  grb2_fold_Delsub$type<-"delSub"    
  grb2_fold_allins$type<-"allins"
  
  ## data normalisation
  # normalise the data to mode of dead peak of the bimodal distribution of variants
  
  #bind all data
  test<-rbind(grb2_fold_synonymous,
              grb2_fold_allins,
              grb2_fold_singles,
              grb2_fold_CX,
              grb2_fold_CXX,
              grb2_fold_CNN,
              grb2_fold_insAA,
              grb2_fold_insCC,
              grb2_fold_CXXX,
              grb2_fold_CNNN,
              grb2_fold_singleDEL,
              grb2_fold_doubleDEL, 
              grb2_fold_tripleDEL,
              fill=TRUE)
  
  ##find mode of the lower peak in the bimodal distribution --> referred to as STOPs
  STOPs<-min(mlv(na.omit(test$fitness), method="naive"))

  ## substract the mode from all variants. 
  grb2_fold_CX$norm_fitness<-grb2_fold_CX$fitness-STOPs
  grb2_fold_CXX$norm_fitness<-grb2_fold_CXX$fitness-STOPs
  grb2_fold_CXXX$norm_fitness<-grb2_fold_CXXX$fitness-STOPs
  grb2_fold_singleDEL$norm_fitness<-grb2_fold_singleDEL$fitness-STOPs
  grb2_fold_doubleDEL$norm_fitness<-grb2_fold_doubleDEL$fitness-STOPs
  grb2_fold_tripleDEL$norm_fitness<-grb2_fold_tripleDEL$fitness-STOPs
  grb2_fold_Delsub$norm_fitness<-grb2_fold_Delsub$fitness-STOPs
  grb2_fold_synonymous$norm_fitness<-grb2_fold_synonymous$fitness-STOPs
  grb2_fold_insAA$norm_fitness<-grb2_fold_insAA$fitness-STOPs
  grb2_fold_insCC$norm_fitness<-grb2_fold_insAA$fitness-STOPs
  grb2_fold_CNN$norm_fitness<-grb2_fold_CNN$fitness-STOPs
  grb2_fold_CNNN$norm_fitness<-grb2_fold_CNNN$fitness-STOPs
  grb2_fold_allins$norm_fitness<-grb2_fold_allins$fitness-STOPs
  grb2_fold_singles$norm_fitness<-grb2_fold_singles$fitness-STOPs
  
  ## find weighted mean of the synonymous variants
  syns<-sum(grb2_fold_synonymous$norm_fitness/(grb2_fold_synonymous$sigma^2))/sum(1/(grb2_fold_synonymous$sigma^2))
  
  ##normalise by dividing variants with the weighted mean of the synonymous variants
  grb2_fold_synonymous$scaled_fitness<-grb2_fold_synonymous$norm_fitness/syns
  grb2_fold_CX$scaled_fitness<-grb2_fold_CX$norm_fitness/syns
  grb2_fold_CXX$scaled_fitness<-grb2_fold_CXX$norm_fitness/syns
  grb2_fold_CXXX$scaled_fitness<-grb2_fold_CXXX$norm_fitness/syns
  grb2_fold_singleDEL$scaled_fitness<-grb2_fold_singleDEL$norm_fitness/syns
  grb2_fold_doubleDEL$scaled_fitness<-grb2_fold_doubleDEL$norm_fitness/syns
  grb2_fold_tripleDEL$scaled_fitness<-grb2_fold_tripleDEL$norm_fitness/syns
  grb2_fold_Delsub$scaled_fitness<-grb2_fold_Delsub$norm_fitness/syns
  grb2_fold_insAA$scaled_fitness<-grb2_fold_insAA$norm_fitness/syns
  grb2_fold_insCC$scaled_fitness<-grb2_fold_insCC$norm_fitness/syns
  grb2_fold_CNN$scaled_fitness<-grb2_fold_CNN$norm_fitness/syns
  grb2_fold_CNNN$scaled_fitness<-grb2_fold_CNNN$norm_fitness/syns
  grb2_fold_allins$scaled_fitness<-grb2_fold_allins$norm_fitness/syns
  grb2_fold_singles$scaled_fitness<-grb2_fold_singles$norm_fitness/syns
  
  ## scale the error by dividing by the same factor: dividing by the weighted mean of the synonmous fitness
  grb2_fold_synonymous$scaled_sigma<-grb2_fold_synonymous$sigma/(syns)
  grb2_fold_CX$scaled_sigma<-grb2_fold_CX$sigma/(syns)
  grb2_fold_CXX$scaled_sigma<-grb2_fold_CXX$sigma/(syns)
  grb2_fold_CXXX$scaled_sigma<-grb2_fold_CXXX$sigma/(syns)
  grb2_fold_singleDEL$scaled_sigma<-grb2_fold_singleDEL$sigma/(syns)
  grb2_fold_doubleDEL$scaled_sigma<-grb2_fold_doubleDEL$sigma/(syns)
  grb2_fold_tripleDEL$scaled_sigma<-grb2_fold_tripleDEL$sigma/(syns)
  grb2_fold_Delsub$scaled_sigma<-grb2_fold_Delsub$sigma/(syns)
  grb2_fold_insAA$scaled_sigma<-grb2_fold_insAA$sigma/(syns)
  grb2_fold_insCC$scaled_sigma<-grb2_fold_insCC$sigma/(syns)
  grb2_fold_CNN$scaled_sigma<-grb2_fold_CNN$sigma/(syns)
  grb2_fold_CNNN$scaled_sigma<-grb2_fold_CNNN$sigma/(syns)
  grb2_fold_allins$scaled_sigma<-grb2_fold_allins$sigma/(syns)
  grb2_fold_singles$scaled_sigma<-grb2_fold_singles$sigma/(syns)
  
  ##fix mutation positions that are called wrong manually. 
  # this is just an artifact of the function used to call the positions of difference and differing aa sequence
  grb2_fold_doubleDEL[grb2_fold_doubleDEL$diff == "c(20, 22)", ]$Pos<-21
  grb2_fold_doubleDEL[grb2_fold_doubleDEL$diff == "c(35, 37)", ]$Pos<-36
  
  grb2_fold_tripleDEL[grb2_fold_tripleDEL$diff == "c(7, 8, 10)", ]$Pos<-8
  grb2_fold_tripleDEL[grb2_fold_tripleDEL$diff == "c(7, 8, 11)", ]$Pos<-9
  grb2_fold_tripleDEL[grb2_fold_tripleDEL$diff == "c(20, 22, 23)", ]$Pos<-21
  grb2_fold_tripleDEL[grb2_fold_tripleDEL$diff == "c(35, 37, 38)", ]$Pos<-36
  
  ## add missing variants to the df
  # these are variants that have the same aa_seq although the mutation is technically in different positons. 
  #singleDEL:20=21,35=36
  grb2_fold_singleDEL_missing<-grb2_fold_singleDEL[c(40,18),]
  grb2_fold_singleDEL_missing$Pos<-c(21,36)
  grb2_fold_singleDEL<-rbind(grb2_fold_singleDEL,
                             grb2_fold_singleDEL_missing)
  
  #doubleDEL: 7=8, 7=9
  grb2_fold_doubleDEL_missing<-grb2_fold_doubleDEL[c(3,3),]
  grb2_fold_doubleDEL_missing$Pos<-c(8,9)
  grb2_fold_doubleDEL<-rbind(grb2_fold_doubleDEL,
                             grb2_fold_doubleDEL_missing)
  
  #tripleDEL:13=14, 15=16, 29=30, 42=43
  grb2_fold_tripleDEL_missing<-grb2_fold_tripleDEL[c(8, 41, 16, 24),]
  grb2_fold_tripleDEL_missing$Pos<-c(14, 16, 30,43)
  grb2_fold_tripleDEL<-rbind(grb2_fold_tripleDEL,
                             grb2_fold_tripleDEL_missing)
  #CX 20=21, 35=36
  grb2_fold_CX_missing<-grb2_fold_CX[c(11,35),]
  grb2_fold_CX_missing$Pos<-c(21,36)
  grb2_fold_CX<-rbind(grb2_fold_CX,
                      grb2_fold_CX_missing)
  
  #### test if the change in fitness is signifcant
  ## we need to test for significance for values > and < than 1.  
  ## We calculate a z-stat and then do a 2-tailed test
  # we do mulitple testing correction using bonferroni
  
  # CX
  zstats=(grb2_fold_CX$scaled_fitness-1)/grb2_fold_CX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_fold_CX$zstat=zstats
  grb2_fold_CX$pvalue=pvals_man
  grb2_fold_CX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_fold_CX$fdrsig=grb2_fold_CX$FDR<=.05
  grb2_fold_CX<-grb2_fold_CX %>% mutate(significant=
                                          ifelse(fdrsig == "TRUE", "*",NA))
  
  
  ##CXX
  zstats=(grb2_fold_CXX$scaled_fitness-1)/grb2_fold_CXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_fold_CXX$zstat=zstats
  grb2_fold_CXX$pvalue=pvals_man
  grb2_fold_CXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_fold_CXX$fdrsig=grb2_fold_CXX$FDR<=.05
  grb2_fold_CXX<-grb2_fold_CXX %>% mutate(significant=
                                            ifelse(fdrsig == "TRUE", "*",NA))
  ##CXXX
  zstats=(grb2_fold_CXXX$scaled_fitness-1)/grb2_fold_CXXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_fold_CXXX$zstat=zstats
  grb2_fold_CXXX$pvalue=pvals_man
  grb2_fold_CXXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_fold_CXXX$fdrsig=grb2_fold_CXXX$FDR<=.05
  grb2_fold_CXXX<-grb2_fold_CXXX %>% mutate(significant=
                                              ifelse(fdrsig == "TRUE", "*",NA))
  #singleDEL
  zstats=(grb2_fold_singleDEL$scaled_fitness-1)/grb2_fold_singleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_fold_singleDEL$zstat=zstats
  grb2_fold_singleDEL$pvalue=pvals_man
  grb2_fold_singleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_fold_singleDEL$fdrsig=grb2_fold_singleDEL$FDR<=.05
  grb2_fold_singleDEL<-grb2_fold_singleDEL %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #doubleDEL
  zstats=(grb2_fold_doubleDEL$scaled_fitness-1)/grb2_fold_doubleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_fold_doubleDEL$zstat=zstats
  grb2_fold_doubleDEL$pvalue=pvals_man
  grb2_fold_doubleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_fold_doubleDEL$fdrsig=grb2_fold_doubleDEL$FDR<=.05
  grb2_fold_doubleDEL<-grb2_fold_doubleDEL %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #tripleDEL
  zstats=(grb2_fold_tripleDEL$scaled_fitness-1)/grb2_fold_tripleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_fold_tripleDEL$zstat=zstats
  grb2_fold_tripleDEL$pvalue=pvals_man
  grb2_fold_tripleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_fold_tripleDEL$fdrsig=grb2_fold_tripleDEL$FDR<=.05
  grb2_fold_tripleDEL<-grb2_fold_tripleDEL %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #DelSub
  zstats=(grb2_fold_Delsub$scaled_fitness-1)/grb2_fold_Delsub$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_fold_Delsub$zstat=zstats
  grb2_fold_Delsub$pvalue=pvals_man
  grb2_fold_Delsub[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_fold_Delsub$fdrsig=grb2_fold_Delsub$FDR<=.05
  grb2_fold_Delsub<-grb2_fold_Delsub %>% mutate(significant=
                                                  ifelse(fdrsig == "TRUE", "*",NA))
  
  grb2_fold_Delsub<-grb2_fold_Delsub[grb2_fold_Delsub$STOP==FALSE,]
  
  
  ##some of the DelSubs can cause 2 different substitutions at the same positions. We need to seperate them out into 2 dataframes. 
  grb2_fold_Delsub_2<-grb2_fold_Delsub[duplicated(grb2_fold_Delsub$Pos),]
  grb2_fold_Delsub<-grb2_fold_Delsub[!duplicated(grb2_fold_Delsub$Pos),]
  grb2_fold_Delsub_2$type<-"delSub_1"
  grb2_fold_Delsub$type<-"delSub_2"
  
  #### CNN
  zstats=(grb2_fold_CNN$scaled_fitness-1)/grb2_fold_CNN$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_fold_CNN$zstat=zstats
  grb2_fold_CNN$pvalue=pvals_man
  grb2_fold_CNN[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_fold_CNN$fdrsig=grb2_fold_CNN$FDR<=.05
  grb2_fold_CNN<-grb2_fold_CNN %>% mutate(significant=
                                            ifelse(fdrsig == "TRUE", "*",NA))
  
  grb2_fold_CNN$Pos<-str_first_number(grb2_fold_CNN$diff)
  
  #### insAA
  zstats=(grb2_fold_insAA$scaled_fitness-1)/grb2_fold_insAA$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_fold_insAA$zstat=zstats
  grb2_fold_insAA$pvalue=pvals_man
  grb2_fold_insAA[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_fold_insAA$fdrsig=grb2_fold_insAA$FDR<=.05
  grb2_fold_insAA<-grb2_fold_insAA %>% mutate(significant=
                                                ifelse(fdrsig == "TRUE", "*",NA))
  
  grb2_fold_insAA$Pos<-str_first_number(grb2_fold_insAA$diff)
  grb2_fold_insAA$type<-"doubleAA"
  
  ### insCC
  zstats=(grb2_fold_insCC$scaled_fitness-1)/grb2_fold_insCC$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_fold_insCC$zstat=zstats
  grb2_fold_insCC$pvalue=pvals_man
  grb2_fold_insCC[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_fold_insCC$fdrsig=grb2_fold_insCC$FDR<=.05
  grb2_fold_insCC<-grb2_fold_insCC %>% mutate(significant=
                                                ifelse(fdrsig == "TRUE", "*",NA))
  
  grb2_fold_insCC$Pos<-str_first_number(grb2_fold_insCC$diff)
  grb2_fold_insCC$type<-"doubleCC"
  
  ### CNNN
  zstats=(grb2_fold_CNNN$scaled_fitness-1)/grb2_fold_CNNN$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_fold_CNNN$zstat=zstats
  grb2_fold_CNNN$pvalue=pvals_man
  grb2_fold_CNNN[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_fold_CNNN$fdrsig=grb2_fold_CNNN$FDR<=.05
  grb2_fold_CNNN<-grb2_fold_CNNN %>% mutate(significant=
                                              ifelse(fdrsig == "TRUE", "*",NA))
  
  grb2_fold_CNNN$Pos<-str_first_number(grb2_fold_CNNN$diff)
  
  ### allins
  zstats=(grb2_fold_allins$scaled_fitness-1)/grb2_fold_allins$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_fold_allins$zstat=zstats
  grb2_fold_allins$pvalue=pvals_man
  grb2_fold_allins[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_fold_allins$fdrsig=grb2_fold_allins$FDR<=.05
  grb2_fold_allins<-grb2_fold_allins %>% mutate(significant=
                                                  ifelse(fdrsig == "TRUE", "*",NA))
  ### singles
  zstats=(grb2_fold_singles$scaled_fitness-1)/grb2_fold_singles$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_fold_singles$zstat=zstats
  grb2_fold_singles$pvalue=pvals_man
  grb2_fold_singles[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_fold_singles$fdrsig=grb2_fold_singles$FDR<=.05
  grb2_fold_singles<-grb2_fold_singles %>% mutate(significant=
                                                    ifelse(fdrsig == "TRUE", "*",NA))
  grb2_fold_singles$type<-"substitution"

  ## bind all variants. 
  grb2_fold_allvariants<-rbind(grb2_fold_synonymous[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
                               grb2_fold_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               grb2_fold_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               grb2_fold_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               grb2_fold_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               grb2_fold_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")], 
                               grb2_fold_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               grb2_fold_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               grb2_fold_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               grb2_fold_singles[grb2_fold_singles$STOP=="FALSE",c("aa_seq", "Pos", "Mut","type", "scaled_fitness", "scaled_sigma", "significant")],
                               fill=TRUE)
  
  
  return(list(grb2_fold_allvariants = grb2_fold_allvariants, 
              grb2_fold_CX = grb2_fold_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_fold_CXX = grb2_fold_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_fold_CXXX = grb2_fold_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_fold_CNN = grb2_fold_CNN[,c("aa_seq", "Pos","insID","lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_fold_CNNN = grb2_fold_CNNN[,c("aa_seq", "Pos", "insID","lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_fold_insCC = grb2_fold_insCC[,c("aa_seq", "Pos", "insID","lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_fold_insAA = grb2_fold_insAA[,c("aa_seq", "Pos", "insID","lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_fold_allins = grb2_fold_allins[,c("aa_seq", "Pos", "insID","lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_fold_synonymous = grb2_fold_synonymous[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
              grb2_fold_singles = grb2_fold_singles[,c("aa_seq", "Pos", "Mut", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_fold_Delsub = grb2_fold_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_fold_Delsub_2 = grb2_fold_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_fold_singleDEL = grb2_fold_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_fold_doubleDEL = grb2_fold_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_fold_tripleDEL = grb2_fold_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")]))
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
