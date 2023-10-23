normalise_P61024_PF01111_data <- function() {
  
  ##isolate only the programmed lenghts and split into different dfs depending on indel type
  P61024_PF01111_variants$lenght=str_length(P61024_PF01111_variants$aa_seq)
  
  P61024_PF01111_singleDEL <- P61024_PF01111_variants[P61024_PF01111_variants$lenght == "66",]
  P61024_PF01111_doubleDEL <- P61024_PF01111_variants[P61024_PF01111_variants$lenght == "65",]
  P61024_PF01111_tripleDEL <- P61024_PF01111_variants[P61024_PF01111_variants$lenght == "64",]
  P61024_PF01111_CX <- P61024_PF01111_variants[P61024_PF01111_variants$lenght == "68",]
  P61024_PF01111_CXX <- P61024_PF01111_variants[P61024_PF01111_variants$lenght == "69",]
  P61024_PF01111_CXXX <- P61024_PF01111_variants[P61024_PF01111_variants$lenght == "70",]
  
  P61024_PF01111_singleDEL<-P61024_PF01111_singleDEL[P61024_PF01111_singleDEL$characters<3,]
  
  ##Delsubs
  P61024_PF01111_Delsub<-P61024_PF01111_variants[P61024_PF01111_variants$lenght == "66",]
  P61024_PF01111_Delsub<-P61024_PF01111_Delsub[P61024_PF01111_Delsub$characters>2,]
  
  ##extract position of mutation, add Pos  and length column
  P61024_PF01111_syns$diff<-as.character(P61024_PF01111_syns$diff)
  P61024_PF01111_syns$diff<-str_first_number(P61024_PF01111_syns$diff)
  P61024_PF01111_syns$diff<-as.numeric(P61024_PF01111_syns$diff)
  P61024_PF01111_syns$Pos<-as.numeric(P61024_PF01111_syns$diff)
  P61024_PF01111_syns$lenght<-str_length(P61024_PF01111_syns$aa_seq)
  
  ##extract positions as numeric into the Pos column for all df
  P61024_PF01111_singleDEL$Pos<-as.numeric(P61024_PF01111_singleDEL$diff)
  
  P61024_PF01111_doubleDEL$diff<-as.character(P61024_PF01111_doubleDEL$diff)
  P61024_PF01111_doubleDEL$Pos<-str_first_number(P61024_PF01111_doubleDEL$diff)
  P61024_PF01111_doubleDEL$Pos<-as.numeric(P61024_PF01111_doubleDEL$Pos)
  
  P61024_PF01111_tripleDEL$diff<-as.character(P61024_PF01111_tripleDEL$diff)
  P61024_PF01111_tripleDEL$Pos<-str_first_number(P61024_PF01111_tripleDEL$diff)
  P61024_PF01111_tripleDEL$Pos<-as.numeric(P61024_PF01111_tripleDEL$Pos)
  
  P61024_PF01111_CX$Pos<-as.numeric(P61024_PF01111_CX$diff)
  
  P61024_PF01111_CXX$diff<-as.character(P61024_PF01111_CXX$diff)
  P61024_PF01111_CXX$Pos<-str_first_number(P61024_PF01111_CXX$diff)
  P61024_PF01111_CXX$Pos<-as.numeric(P61024_PF01111_CXX$Pos)
  
  P61024_PF01111_CXXX$diff<-as.character(P61024_PF01111_CXXX$diff)
  P61024_PF01111_CXXX$Pos<-str_first_number(P61024_PF01111_CXXX$diff)
  P61024_PF01111_CXXX$Pos<-as.numeric(P61024_PF01111_CXXX$Pos)
  
  P61024_PF01111_Delsub$diff<-as.character(P61024_PF01111_Delsub$diff)
  P61024_PF01111_Delsub$Pos<-str_first_number(P61024_PF01111_Delsub$diff)
  
  ### find the aa position of the synonymous mutation
  P61024_PF01111_syns$Pos <- as.character(floor(P61024_PF01111_syns$diff / 3))
  
  
  ## add mutation type
  P61024_PF01111_singleDEL$type<-"singleDEL"
  P61024_PF01111_doubleDEL$type<-"doubleDEL"
  P61024_PF01111_tripleDEL$type<-"tripleDEL"
  P61024_PF01111_CX$type<-"singleINS"
  P61024_PF01111_CXX$type<-"doubleINS"
  P61024_PF01111_CXXX$type<-"tripleINS"
  P61024_PF01111_Delsub$type<-"delSub"    
  P61024_PF01111_syns$type<-"synonymous"
  
  ##fix mutation positions that are called wrong manually. 
  # this is just an artifact of the function used to call the positions of difference and differing aa sequence
  P61024_PF01111_CXX[P61024_PF01111_CXX$diff == "c(8, 10)",]$Pos<-9
  P61024_PF01111_CXX[P61024_PF01111_CXX$diff == "c(10, 12)",]$Pos<-11
  P61024_PF01111_CXX[P61024_PF01111_CXX$diff == "c(44, 46)",]$Pos<-45
  P61024_PF01111_CXX[P61024_PF01111_CXX$diff == "c(62, 64)",]$Pos<-63
  P61024_PF01111_CXX[P61024_PF01111_CXX$diff == "c(65, 67)",]$Pos<-66
  
  P61024_PF01111_CXXX[P61024_PF01111_CXXX$diff == "c(2, 4, 5)",]$Pos<-3
  P61024_PF01111_CXXX[P61024_PF01111_CXXX$diff == "c(8, 10, 11)",]$Pos<-9
  P61024_PF01111_CXXX[P61024_PF01111_CXXX$diff == "c(10, 12, 14)",]$Pos<-12
  P61024_PF01111_CXXX[P61024_PF01111_CXXX$diff == "c(34, 35, 38)",]$Pos<-35
  P61024_PF01111_CXXX[P61024_PF01111_CXXX$diff == "c(44, 46, 48)",]$Pos<-46
  P61024_PF01111_CXXX[P61024_PF01111_CXXX$diff == "c(56, 57, 59)",]$Pos<-57
  P61024_PF01111_CXXX[P61024_PF01111_CXXX$diff == "c(56, 57, 60)",]$Pos<-58
  
  P61024_PF01111_doubleDEL[P61024_PF01111_doubleDEL$diff == "c(8, 10)",]$Pos<-9
  P61024_PF01111_doubleDEL[P61024_PF01111_doubleDEL$diff == "c(10, 12)",]$Pos<-12
  P61024_PF01111_doubleDEL[P61024_PF01111_doubleDEL$diff == "c(44, 46)",]$Pos<-46
  P61024_PF01111_doubleDEL[P61024_PF01111_doubleDEL$diff == "c(62, 64)",]$Pos<-63
  P61024_PF01111_doubleDEL[P61024_PF01111_doubleDEL$diff == "c(65, 67)",]$Pos<-66
  
  P61024_PF01111_tripleDEL[P61024_PF01111_tripleDEL$diff == "c(2, 4, 5)",]$Pos<-3
  P61024_PF01111_tripleDEL[P61024_PF01111_tripleDEL$diff == "c(10, 12, 14)",]$Pos<-12
  P61024_PF01111_tripleDEL[P61024_PF01111_tripleDEL$diff == "c(34, 35, 37)",]$Pos<-35
  P61024_PF01111_tripleDEL[P61024_PF01111_tripleDEL$diff == "c(34, 35, 38)",]$Pos<-36
  P61024_PF01111_tripleDEL[P61024_PF01111_tripleDEL$diff == "c(44, 46, 48)",]$Pos<-45
  P61024_PF01111_tripleDEL[P61024_PF01111_tripleDEL$diff == "c(56, 57, 59)",]$Pos<-57
  P61024_PF01111_tripleDEL[P61024_PF01111_tripleDEL$diff == "c(56, 57, 60)",]$Pos<-58
  P61024_PF01111_tripleDEL[P61024_PF01111_tripleDEL$diff == "c(62, 64, 65)",]$Pos<-63
  
  # add missing variants to the df
  # these are variants that have the same aa_seq although the mutation is technically in different positons. 
  #singleDEL: 2=3, 8=9, 10=11, 44=45, 62=63, 65=66
  P61024_PF01111_singleDEL_missing<-P61024_PF01111_singleDEL[c(50, 4, 48, 28, 27),]
  P61024_PF01111_singleDEL_missing$Pos<-c(3, 9, 11, 63, 66)
  P61024_PF01111_singleDEL<-rbind(P61024_PF01111_singleDEL,
                                  P61024_PF01111_singleDEL_missing)
  #doubleDEL: 11=12, 34=35, 34=36, 45=46, 56=57, 56=58
  P61024_PF01111_doubleDEL_missing<-P61024_PF01111_doubleDEL[c(6, 42, 42, 21,25,25),]
  P61024_PF01111_doubleDEL_missing$Pos<-c(11,35,36,45,57,58)
  P61024_PF01111_doubleDEL<-rbind(P61024_PF01111_doubleDEL,
                                  P61024_PF01111_doubleDEL_missing)
  #tripleDEL:5=6, 10=11, 44=45
  P61024_PF01111_tripleDEL_missing<-P61024_PF01111_tripleDEL[c(57, 56, 41),]
  P61024_PF01111_tripleDEL_missing$Pos<-c(6, 11, 45)
  P61024_PF01111_tripleDEL<-rbind(P61024_PF01111_tripleDEL,
                                  P61024_PF01111_tripleDEL_missing)
  #CX 2=3, 8=9, 10=11, 44=45, 62=63, 65=66
  P61024_PF01111_CX_missing<-P61024_PF01111_CX[c(2, 48, 3, 15, 23, 24),]
  P61024_PF01111_CX_missing$Pos<-c(3,9,11,45,63,66)
  P61024_PF01111_CX<-rbind(P61024_PF01111_CX,
                           P61024_PF01111_CX_missing)
  #CXX 11=12, 34=35, 34=36, 45=46, 56=57, 56=58
  P61024_PF01111_CXX_missing<-P61024_PF01111_CXX[c(45, 13,13, 34, 30,30),]
  P61024_PF01111_CXX_missing$Pos<-c(12,35,36,46,57,58)
  P61024_PF01111_CXX<-rbind(P61024_PF01111_CXX,
                            P61024_PF01111_CXX_missing)
  #CXXX 5=6, 10=11, 44=45
  P61024_PF01111_CXXX_missing<-P61024_PF01111_CXXX[c(1, 13),]
  P61024_PF01111_CXXX_missing$Pos<-c(6,45)
  P61024_PF01111_CXXX<-rbind(P61024_PF01111_CXXX,
                             P61024_PF01111_CXXX_missing)
  
  ## data normalisation
  # normalise the data to mode of dead peak of the bimodal distribution of variants
  #bind all data
  test<-rbind(P61024_PF01111_syns,
              P61024_PF01111_CX,
              P61024_PF01111_CXX,
              P61024_PF01111_CXXX,
              P61024_PF01111_singleDEL,
              P61024_PF01111_doubleDEL, 
              P61024_PF01111_tripleDEL,
              fill=TRUE)
  
  ##find mode of the lower peak in the bimodal distribution --> referred to as STOPs
  STOPs<-mlv(na.omit(test$fitness), method="naive")

  ## substract the mode from all variants. 
  P61024_PF01111_CX$norm_fitness<-P61024_PF01111_CX$fitness-STOPs
  P61024_PF01111_CXX$norm_fitness<-P61024_PF01111_CXX$fitness-STOPs
  P61024_PF01111_CXXX$norm_fitness<-P61024_PF01111_CXXX$fitness-STOPs
  P61024_PF01111_singleDEL$norm_fitness<-P61024_PF01111_singleDEL$fitness-STOPs
  P61024_PF01111_doubleDEL$norm_fitness<-P61024_PF01111_doubleDEL$fitness-STOPs
  P61024_PF01111_tripleDEL$norm_fitness<-P61024_PF01111_tripleDEL$fitness-STOPs
  P61024_PF01111_Delsub$norm_fitness<-P61024_PF01111_Delsub$fitness-STOPs
  P61024_PF01111_syns$norm_fitness<-P61024_PF01111_syns$fitness-STOPs
  
  ## find weighted mean of the synonymous variants
  syns<-sum(P61024_PF01111_syns$norm_fitness/(P61024_PF01111_syns$sigma^2))/sum(1/(P61024_PF01111_syns$sigma^2))
  
  ##normalise by dividing variants with the weighted mean of the synonymous variants
  P61024_PF01111_syns$scaled_fitness<-P61024_PF01111_syns$norm_fitness/syns
  P61024_PF01111_CX$scaled_fitness<-P61024_PF01111_CX$norm_fitness/syns
  P61024_PF01111_CXX$scaled_fitness<-P61024_PF01111_CXX$norm_fitness/syns
  P61024_PF01111_CXXX$scaled_fitness<-P61024_PF01111_CXXX$norm_fitness/syns
  P61024_PF01111_singleDEL$scaled_fitness<-P61024_PF01111_singleDEL$norm_fitness/syns
  P61024_PF01111_doubleDEL$scaled_fitness<-P61024_PF01111_doubleDEL$norm_fitness/syns
  P61024_PF01111_tripleDEL$scaled_fitness<-P61024_PF01111_tripleDEL$norm_fitness/syns
  P61024_PF01111_Delsub$scaled_fitness<-P61024_PF01111_Delsub$norm_fitness/syns
  
  
  ## scale the error by dividing by the same factor: dividing by the weighted mean of the synonmous fitness
  P61024_PF01111_syns$scaled_sigma<-P61024_PF01111_syns$sigma/(syns)
  P61024_PF01111_CX$scaled_sigma<-P61024_PF01111_CX$sigma/(syns)
  P61024_PF01111_CXX$scaled_sigma<-P61024_PF01111_CXX$sigma/(syns)
  P61024_PF01111_CXXX$scaled_sigma<-P61024_PF01111_CXXX$sigma/(syns)
  P61024_PF01111_singleDEL$scaled_sigma<-P61024_PF01111_singleDEL$sigma/(syns)
  P61024_PF01111_doubleDEL$scaled_sigma<-P61024_PF01111_doubleDEL$sigma/(syns)
  P61024_PF01111_tripleDEL$scaled_sigma<-P61024_PF01111_tripleDEL$sigma/(syns)
  P61024_PF01111_Delsub$scaled_sigma<-P61024_PF01111_Delsub$sigma/(syns)

  #### test if the change in fitness is signifcant
  ## we need to test for significance for values > and < than 1.  
  ## We calculate a z-stat and then do a 2-tailed test
  # we do mulitple testing correction using bonferroni
  
  # CX
  zstats=(P61024_PF01111_CX$scaled_fitness-1)/P61024_PF01111_CX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P61024_PF01111_CX$zstat=zstats
  P61024_PF01111_CX$pvalue=pvals_man
  P61024_PF01111_CX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P61024_PF01111_CX$fdrsig=P61024_PF01111_CX$FDR<=.05
  P61024_PF01111_CX<-P61024_PF01111_CX %>% mutate(significant=
                                                    ifelse(fdrsig == "TRUE", "*",NA))
  
  
  ##CXX
  zstats=(P61024_PF01111_CXX$scaled_fitness-1)/P61024_PF01111_CXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P61024_PF01111_CXX$zstat=zstats
  P61024_PF01111_CXX$pvalue=pvals_man
  P61024_PF01111_CXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P61024_PF01111_CXX$fdrsig=P61024_PF01111_CXX$FDR<=.05
  P61024_PF01111_CXX<-P61024_PF01111_CXX %>% mutate(significant=
                                                      ifelse(fdrsig == "TRUE", "*",NA))
  ##CXXX
  zstats=(P61024_PF01111_CXXX$scaled_fitness-1)/P61024_PF01111_CXXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P61024_PF01111_CXXX$zstat=zstats
  P61024_PF01111_CXXX$pvalue=pvals_man
  P61024_PF01111_CXXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P61024_PF01111_CXXX$fdrsig=P61024_PF01111_CXXX$FDR<=.05
  P61024_PF01111_CXXX<-P61024_PF01111_CXXX %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #singleDEL
  zstats=(P61024_PF01111_singleDEL$scaled_fitness-1)/P61024_PF01111_singleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P61024_PF01111_singleDEL$zstat=zstats
  P61024_PF01111_singleDEL$pvalue=pvals_man
  P61024_PF01111_singleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P61024_PF01111_singleDEL$fdrsig=P61024_PF01111_singleDEL$FDR<=.05
  P61024_PF01111_singleDEL<-P61024_PF01111_singleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #doubleDEL
  zstats=(P61024_PF01111_doubleDEL$scaled_fitness-1)/P61024_PF01111_doubleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P61024_PF01111_doubleDEL$zstat=zstats
  P61024_PF01111_doubleDEL$pvalue=pvals_man
  P61024_PF01111_doubleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P61024_PF01111_doubleDEL$fdrsig=P61024_PF01111_doubleDEL$FDR<=.05
  P61024_PF01111_doubleDEL<-P61024_PF01111_doubleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #tripleDEL
  zstats=(P61024_PF01111_tripleDEL$scaled_fitness-1)/P61024_PF01111_tripleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P61024_PF01111_tripleDEL$zstat=zstats
  P61024_PF01111_tripleDEL$pvalue=pvals_man
  P61024_PF01111_tripleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P61024_PF01111_tripleDEL$fdrsig=P61024_PF01111_tripleDEL$FDR<=.05
  P61024_PF01111_tripleDEL<-P61024_PF01111_tripleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #DelSub
  zstats=(P61024_PF01111_Delsub$scaled_fitness-1)/P61024_PF01111_Delsub$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P61024_PF01111_Delsub$zstat=zstats
  P61024_PF01111_Delsub$pvalue=pvals_man
  P61024_PF01111_Delsub[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P61024_PF01111_Delsub$fdrsig=P61024_PF01111_Delsub$FDR<=.05
  P61024_PF01111_Delsub<-P61024_PF01111_Delsub %>% mutate(significant=
                                                            ifelse(fdrsig == "TRUE", "*",NA))
  
  ##some of the DelSubs can cause 2 different substitutions at the same positions. We need to seperate them out into 2 dataframes. 
  P61024_PF01111_Delsub_2<-P61024_PF01111_Delsub[duplicated(P61024_PF01111_Delsub$Pos),]
  P61024_PF01111_Delsub<-P61024_PF01111_Delsub[!duplicated(P61024_PF01111_Delsub$Pos),]
  P61024_PF01111_Delsub_2$type<-"delSub_1"
  P61024_PF01111_Delsub$type<-"delSub_2"
  P61024_PF01111_Delsub_2$Pos<-as.numeric(P61024_PF01111_Delsub_2$Pos)
  P61024_PF01111_Delsub$Pos<-as.numeric(P61024_PF01111_Delsub$Pos)
  
  P61024_PF01111_Delsub<-P61024_PF01111_Delsub[P61024_PF01111_Delsub$STOP==FALSE,]
  
  
  ##add missing variants with 0 in output:
  P61024_PF01111_missing$lenght=nchar(P61024_PF01111_missing$aa_seq)
  ## remove columns
  P61024_PF01111_missing<-P61024_PF01111_missing[,c(2,21, 23:24)]
  ## isolate the missing variants by lenght 
  P61024_PF01111_singleDEL_missing <- P61024_PF01111_missing[P61024_PF01111_missing$lenght == "66",]
  P61024_PF01111_doubleDEL_missing <- P61024_PF01111_missing[P61024_PF01111_missing$lenght == "65",]
  P61024_PF01111_tripleDEL_missing <- P61024_PF01111_missing[P61024_PF01111_missing$lenght == "64",]
  P61024_PF01111_CX_missing <- P61024_PF01111_missing[P61024_PF01111_missing$lenght == "68",]
  P61024_PF01111_CXX_missing <- P61024_PF01111_missing[P61024_PF01111_missing$lenght == "69",]
  P61024_PF01111_CXXX_missing <- P61024_PF01111_missing[P61024_PF01111_missing$lenght == "70",]
  
  # singleDEL: remove mixed variants (subs and indel, not programmed)
  rows_to_del<-which(nchar(P61024_PF01111_singleDEL_missing$P61024_PF01111)>2)
  P61024_PF01111_singleDEL_missing<-P61024_PF01111_singleDEL_missing[-rows_to_del,]
  P61024_PF01111_singleDEL_missing$Pos<-as.numeric(P61024_PF01111_singleDEL_missing$P61024_PF01111)
  
  P61024_PF01111_singleDEL<-rbind(P61024_PF01111_singleDEL,
                                  P61024_PF01111_singleDEL_missing,
                                  fill=TRUE)
  
  P61024_PF01111_singleDEL$type<-"singleDEL"
  
  # doubleDEL: remove mixed variants (subs and indel, not programmed)
  P61024_PF01111_doubleDEL_missing$Pos<-c(14,10)
  P61024_PF01111_doubleDEL<-rbind(P61024_PF01111_doubleDEL,
                                  P61024_PF01111_doubleDEL_missing,
                                  fill=TRUE)
  
  P61024_PF01111_doubleDEL$type<-"doubleDEL"
  
  # tripleDEL: remove mixed variants (subs and indel, not programmed)
  P61024_PF01111_tripleDEL_missing$Pos<-c(14,55,33)
  P61024_PF01111_tripleDEL<-rbind(P61024_PF01111_tripleDEL,
                                  P61024_PF01111_tripleDEL_missing,
                                  fill=TRUE)
  
  P61024_PF01111_tripleDEL$type<-"tripleDEL"
  
  # singleINS: remove mixed variants (subs and indel, not programmed)
  P61024_PF01111_CX_missing$Pos<-as.numeric(P61024_PF01111_CX_missing$P61024_PF01111)
  P61024_PF01111_CX<-rbind(P61024_PF01111_CX,
                           P61024_PF01111_CX_missing,
                           fill=TRUE)
  
  P61024_PF01111_CX$type<-"singleINS"
  
  # doubleINS: remove mixed variants (subs and indel, not programmed)
  P61024_PF01111_CXX_missing$Pos<-c(18,24,44,41,32,31,16)
  P61024_PF01111_CXX<-rbind(P61024_PF01111_CXX,
                            P61024_PF01111_CXX_missing,
                            fill=TRUE)
  
  P61024_PF01111_CXX$type<-"doubleINS"
  
  # tripleINS: remove mixed variants (subs and indel, not programmed)
  P61024_PF01111_CXXX_missing$Pos<-c(15,24,61,52,51,50,49,26,13)
  P61024_PF01111_CXXX<-rbind(P61024_PF01111_CXXX,
                             P61024_PF01111_CXXX_missing,
                             fill=TRUE)
  
  P61024_PF01111_CXXX$type<-"tripleINS"
  
  ## final df
  P61024_PF01111_allvariants<-rbind(P61024_PF01111_syns[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
                                    P61024_PF01111_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P61024_PF01111_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P61024_PF01111_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P61024_PF01111_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P61024_PF01111_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")], 
                                    P61024_PF01111_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P61024_PF01111_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P61024_PF01111_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    fill=TRUE)
  
  return(list(P61024_PF01111_allvariants = P61024_PF01111_allvariants, 
              P61024_PF01111_syns = P61024_PF01111_syns[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
              P61024_PF01111_CX = P61024_PF01111_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P61024_PF01111_CXX = P61024_PF01111_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P61024_PF01111_CXXX = P61024_PF01111_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P61024_PF01111_singleDEL = P61024_PF01111_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P61024_PF01111_doubleDEL = P61024_PF01111_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P61024_PF01111_tripleDEL = P61024_PF01111_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P61024_PF01111_Delsub = P61024_PF01111_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P61024_PF01111_Delsub_2 = P61024_PF01111_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")]))
  
}







  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
