normalise_P32081_PF00313_data <- function() {
  
  ##isolate only the programmed lenghts and split into different dfs depending on indel type
  P32081_PF00313_variants$lenght=str_length(P32081_PF00313_variants$aa_seq)
  
  P32081_PF00313_singleDEL <- P32081_PF00313_variants[P32081_PF00313_variants$lenght == "62",]
  P32081_PF00313_doubleDEL <- P32081_PF00313_variants[P32081_PF00313_variants$lenght == "61",]
  P32081_PF00313_tripleDEL <- P32081_PF00313_variants[P32081_PF00313_variants$lenght == "60",]
  P32081_PF00313_CX <- P32081_PF00313_variants[P32081_PF00313_variants$lenght == "64",]
  P32081_PF00313_CXX <- P32081_PF00313_variants[P32081_PF00313_variants$lenght == "65",]
  P32081_PF00313_CXXX <- P32081_PF00313_variants[P32081_PF00313_variants$lenght == "66",]
  
  P32081_PF00313_singleDEL<-P32081_PF00313_singleDEL[P32081_PF00313_singleDEL$characters<3,]
  
  ##Delsubs
  P32081_PF00313_Delsub<-P32081_PF00313_variants[P32081_PF00313_variants$lenght == "62",]
  P32081_PF00313_Delsub<-P32081_PF00313_Delsub[P32081_PF00313_Delsub$characters>2,]
  
  ##extract position of mutation, add Pos  and length column
  P32081_PF00313_syns$diff<-as.character(P32081_PF00313_syns$diff)
  P32081_PF00313_syns$diff<-str_first_number(P32081_PF00313_syns$diff)
  P32081_PF00313_syns$diff<-as.numeric(P32081_PF00313_syns$diff)
  P32081_PF00313_syns$Pos<-as.numeric(P32081_PF00313_syns$diff)
  P32081_PF00313_syns$lenght<-str_length(P32081_PF00313_syns$aa_seq)
  
  ##extract positions as numeric into the Pos column for all df
  P32081_PF00313_singleDEL$Pos<-as.numeric(P32081_PF00313_singleDEL$diff)
  
  P32081_PF00313_doubleDEL$diff<-as.character(P32081_PF00313_doubleDEL$diff)
  P32081_PF00313_doubleDEL$Pos<-str_first_number(P32081_PF00313_doubleDEL$diff)
  P32081_PF00313_doubleDEL$Pos<-as.numeric(P32081_PF00313_doubleDEL$Pos)
  
  P32081_PF00313_tripleDEL$diff<-as.character(P32081_PF00313_tripleDEL$diff)
  P32081_PF00313_tripleDEL$Pos<-str_first_number(P32081_PF00313_tripleDEL$diff)
  P32081_PF00313_tripleDEL$Pos<-as.numeric(P32081_PF00313_tripleDEL$Pos)
  
  P32081_PF00313_CX$Pos<-as.numeric(P32081_PF00313_CX$diff)
  
  P32081_PF00313_CXX$diff<-as.character(P32081_PF00313_CXX$diff)
  P32081_PF00313_CXX$Pos<-str_first_number(P32081_PF00313_CXX$diff)
  P32081_PF00313_CXX$Pos<-as.numeric(P32081_PF00313_CXX$Pos)
  
  P32081_PF00313_CXXX$diff<-as.character(P32081_PF00313_CXXX$diff)
  P32081_PF00313_CXXX$Pos<-str_first_number(P32081_PF00313_CXXX$diff)
  P32081_PF00313_CXXX$Pos<-as.numeric(P32081_PF00313_CXXX$Pos)
  
  P32081_PF00313_Delsub$diff<-as.character(P32081_PF00313_Delsub$diff)
  P32081_PF00313_Delsub$Pos<-str_first_number(P32081_PF00313_Delsub$diff)
  
  ### find the aa position of the synonymous mutation
  P32081_PF00313_syns$Pos <- as.character(floor(P32081_PF00313_syns$diff / 3))
  
  
  ## add mutation type
  P32081_PF00313_singleDEL$type<-"singleDEL"
  P32081_PF00313_doubleDEL$type<-"doubleDEL"
  P32081_PF00313_tripleDEL$type<-"tripleDEL"
  P32081_PF00313_CX$type<-"singleINS"
  P32081_PF00313_CXX$type<-"doubleINS"
  P32081_PF00313_CXXX$type<-"tripleINS"
  P32081_PF00313_Delsub$type<-"delSub"    
  P32081_PF00313_syns$type<-"synonymous"
  
  ##fix mutation positions that are called wrong manually. 
  # this is just an artifact of the function used to call the positions of difference and differing aa sequence
  P32081_PF00313_CXX[P32081_PF00313_CXX$diff == "c(22, 24)",]$Pos<-23
  P32081_PF00313_CXX[P32081_PF00313_CXX$diff == "c(40, 42)",]$Pos<-41
  P32081_PF00313_CXX[P32081_PF00313_CXX$diff == "c(58, 60)",]$Pos<-59
  
  P32081_PF00313_CXXX[P32081_PF00313_CXXX$diff == "c(3, 4, 6)",]$Pos<-4
  P32081_PF00313_CXXX[P32081_PF00313_CXXX$diff == "c(12, 13, 15)",]$Pos<-13
  P32081_PF00313_CXXX[P32081_PF00313_CXXX$diff == "c(12, 13, 16)",]$Pos<-14
  P32081_PF00313_CXXX[P32081_PF00313_CXXX$diff == "c(17, 18, 20)",]$Pos<-18
  P32081_PF00313_CXXX[P32081_PF00313_CXXX$diff == "c(22, 24, 25)",]$Pos<-23
  P32081_PF00313_CXXX[P32081_PF00313_CXXX$diff == "c(24, 25, 27)",]$Pos<-25
  P32081_PF00313_CXXX[P32081_PF00313_CXXX$diff == "c(33, 34, 36)",]$Pos<-34
  P32081_PF00313_CXXX[P32081_PF00313_CXXX$diff == "c(40, 42, 43)",]$Pos<-41
  P32081_PF00313_CXXX[P32081_PF00313_CXXX$diff == "c(58, 60, 61)",]$Pos<-59
  
  P32081_PF00313_doubleDEL[P32081_PF00313_doubleDEL$diff == "c(22, 24)",]$Pos<-23
  P32081_PF00313_doubleDEL[P32081_PF00313_doubleDEL$diff == "c(40, 42)",]$Pos<-41
  P32081_PF00313_doubleDEL[P32081_PF00313_doubleDEL$diff == "c(58, 60)",]$Pos<-59
  
  P32081_PF00313_tripleDEL[P32081_PF00313_tripleDEL$diff == "c(3, 4, 6)",]$Pos<-4
  P32081_PF00313_tripleDEL[P32081_PF00313_tripleDEL$diff == "c(12, 13, 15)",]$Pos<-13
  P32081_PF00313_tripleDEL[P32081_PF00313_tripleDEL$diff == "c(12, 13, 16)",]$Pos<-14
  P32081_PF00313_tripleDEL[P32081_PF00313_tripleDEL$diff == "c(17, 18, 20)",]$Pos<-18
  P32081_PF00313_tripleDEL[P32081_PF00313_tripleDEL$diff == "c(22, 24, 25)",]$Pos<-23
  P32081_PF00313_tripleDEL[P32081_PF00313_tripleDEL$diff == "c(24, 25, 27)",]$Pos<-25
  P32081_PF00313_tripleDEL[P32081_PF00313_tripleDEL$diff == "c(33, 34, 36)",]$Pos<-34
  P32081_PF00313_tripleDEL[P32081_PF00313_tripleDEL$diff == "c(40, 42, 43)",]$Pos<-41
  P32081_PF00313_tripleDEL[P32081_PF00313_tripleDEL$diff == "c(58, 60, 61)",]$Pos<-59
  
  ## add missing variants to the df
  # these are variants that have the same aa_seq although the mutation is technically in different positons. 
  #singleDEL: 22=23, 40=41, 58=59
  P32081_PF00313_singleDEL_missing<-P32081_PF00313_singleDEL[c(48, 39, 25),]
  P32081_PF00313_singleDEL_missing$Pos<-c(23,41,59)
  P32081_PF00313_singleDEL<-rbind(P32081_PF00313_singleDEL,
                                  P32081_PF00313_singleDEL_missing)
  #doubleDEL: 3=4, 12=13, 12=14, 17=18, 24=25, 33=34
  P32081_PF00313_doubleDEL_missing<-P32081_PF00313_doubleDEL[c(55, 5, 5, 7, 10, 42),]
  P32081_PF00313_doubleDEL_missing$Pos<-c(4,13,14,18,25,34)
  P32081_PF00313_doubleDEL<-rbind(P32081_PF00313_doubleDEL,
                                  P32081_PF00313_doubleDEL_missing)
  #tripleDEL:25=26, 48=49, 52=53
  P32081_PF00313_tripleDEL_missing<-P32081_PF00313_tripleDEL[c(45, 33, 32),]
  P32081_PF00313_tripleDEL_missing$Pos<-c(26,49,53)
  P32081_PF00313_tripleDEL<-rbind(P32081_PF00313_tripleDEL,
                                  P32081_PF00313_tripleDEL_missing)
  #CX 22=23, 40=41, 58=59
  P32081_PF00313_CX_missing<-P32081_PF00313_CX[c(13, 22, 36),]
  P32081_PF00313_CX_missing$Pos<-c(23,41,59)
  P32081_PF00313_CX<-rbind(P32081_PF00313_CX,
                           P32081_PF00313_CX_missing)
  #CXX 3=4, 12=13, 12=14, 17=18, 24=25, 33=34
  P32081_PF00313_CXX_missing<-P32081_PF00313_CXX[c(2, 52, 52, 50, 47, 15),]
  P32081_PF00313_CXX_missing$Pos<-c(4,13,14,18,25,34)
  P32081_PF00313_CXX<-rbind(P32081_PF00313_CXX,
                            P32081_PF00313_CXX_missing)
  #CXXX 25=26, 48=49, 52=53
  P32081_PF00313_CXXX_missing<-P32081_PF00313_CXXX[c(14, 26, 27),]
  P32081_PF00313_CXXX_missing$Pos<-c(26,49, 53)
  P32081_PF00313_CXXX<-rbind(P32081_PF00313_CXXX,
                             P32081_PF00313_CXXX_missing)
  
  ## data normalisation
  # normalise the data to mode of dead peak of the bimodal distribution of variants
  #bind all data
  test<-rbind(P32081_PF00313_syns,
              P32081_PF00313_CX,
              P32081_PF00313_CXX,
              P32081_PF00313_CXXX,
              P32081_PF00313_singleDEL,
              P32081_PF00313_doubleDEL, 
              P32081_PF00313_tripleDEL,
              fill=TRUE)
  
  ##find mode of the lower peak in the bimodal distribution --> referred to as STOPs
  STOPs<-min(mlv(na.omit(test[test$fitness< -0.2,]$fitness), method="naive"))
  
  ## substract the mode from all variants. 
  P32081_PF00313_CX$norm_fitness<-P32081_PF00313_CX$fitness-STOPs
  P32081_PF00313_CXX$norm_fitness<-P32081_PF00313_CXX$fitness-STOPs
  P32081_PF00313_CXXX$norm_fitness<-P32081_PF00313_CXXX$fitness-STOPs
  P32081_PF00313_singleDEL$norm_fitness<-P32081_PF00313_singleDEL$fitness-STOPs
  P32081_PF00313_doubleDEL$norm_fitness<-P32081_PF00313_doubleDEL$fitness-STOPs
  P32081_PF00313_tripleDEL$norm_fitness<-P32081_PF00313_tripleDEL$fitness-STOPs
  P32081_PF00313_Delsub$norm_fitness<-P32081_PF00313_Delsub$fitness-STOPs
  P32081_PF00313_syns$norm_fitness<-P32081_PF00313_syns$fitness-STOPs
  
  ## find weighted mean of the synonymous variants
  syns<-sum(P32081_PF00313_syns$norm_fitness/(P32081_PF00313_syns$sigma^2))/sum(1/(P32081_PF00313_syns$sigma^2))
  
  ##normalise by dividing variants with the weighted mean of the synonymous variants
  P32081_PF00313_syns$scaled_fitness<-P32081_PF00313_syns$norm_fitness/syns
  P32081_PF00313_CX$scaled_fitness<-P32081_PF00313_CX$norm_fitness/syns
  P32081_PF00313_CXX$scaled_fitness<-P32081_PF00313_CXX$norm_fitness/syns
  P32081_PF00313_CXXX$scaled_fitness<-P32081_PF00313_CXXX$norm_fitness/syns
  P32081_PF00313_singleDEL$scaled_fitness<-P32081_PF00313_singleDEL$norm_fitness/syns
  P32081_PF00313_doubleDEL$scaled_fitness<-P32081_PF00313_doubleDEL$norm_fitness/syns
  P32081_PF00313_tripleDEL$scaled_fitness<-P32081_PF00313_tripleDEL$norm_fitness/syns
  P32081_PF00313_Delsub$scaled_fitness<-P32081_PF00313_Delsub$norm_fitness/syns
  
  
  ## scale the error by dividing by the same factor: dividing by the squared root of the synonmous fitness
  P32081_PF00313_syns$scaled_sigma<-P32081_PF00313_syns$sigma/sqrt(syns)
  P32081_PF00313_CX$scaled_sigma<-P32081_PF00313_CX$sigma/sqrt(syns)
  P32081_PF00313_CXX$scaled_sigma<-P32081_PF00313_CXX$sigma/sqrt(syns)
  P32081_PF00313_CXXX$scaled_sigma<-P32081_PF00313_CXXX$sigma/sqrt(syns)
  P32081_PF00313_singleDEL$scaled_sigma<-P32081_PF00313_singleDEL$sigma/sqrt(syns)
  P32081_PF00313_doubleDEL$scaled_sigma<-P32081_PF00313_doubleDEL$sigma/sqrt(syns)
  P32081_PF00313_tripleDEL$scaled_sigma<-P32081_PF00313_tripleDEL$sigma/sqrt(syns)
  P32081_PF00313_Delsub$scaled_sigma<-P32081_PF00313_Delsub$sigma/sqrt(syns)
  
  #### test if the change in fitness is signifcant
  ## we need to test for significance for values > and < than 1.  
  ## We calculate a z-stat and then do a 2-tailed test
  # we do mulitple testing correction using bonferroni
  
  # CX
  zstats=(P32081_PF00313_CX$scaled_fitness-1)/P32081_PF00313_CX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P32081_PF00313_CX$zstat=zstats
  P32081_PF00313_CX$pvalue=pvals_man
  P32081_PF00313_CX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P32081_PF00313_CX$fdrsig=P32081_PF00313_CX$FDR<=.05
  P32081_PF00313_CX<-P32081_PF00313_CX %>% mutate(significant=
                                                    ifelse(fdrsig == "TRUE", "*",NA))
  
  
  ##CXX
  zstats=(P32081_PF00313_CXX$scaled_fitness-1)/P32081_PF00313_CXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P32081_PF00313_CXX$zstat=zstats
  P32081_PF00313_CXX$pvalue=pvals_man
  P32081_PF00313_CXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P32081_PF00313_CXX$fdrsig=P32081_PF00313_CXX$FDR<=.05
  P32081_PF00313_CXX<-P32081_PF00313_CXX %>% mutate(significant=
                                                      ifelse(fdrsig == "TRUE", "*",NA))
  ##CXXX
  zstats=(P32081_PF00313_CXXX$scaled_fitness-1)/P32081_PF00313_CXXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P32081_PF00313_CXXX$zstat=zstats
  P32081_PF00313_CXXX$pvalue=pvals_man
  P32081_PF00313_CXXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P32081_PF00313_CXXX$fdrsig=P32081_PF00313_CXXX$FDR<=.05
  P32081_PF00313_CXXX<-P32081_PF00313_CXXX %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #singleDEL
  zstats=(P32081_PF00313_singleDEL$scaled_fitness-1)/P32081_PF00313_singleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P32081_PF00313_singleDEL$zstat=zstats
  P32081_PF00313_singleDEL$pvalue=pvals_man
  P32081_PF00313_singleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P32081_PF00313_singleDEL$fdrsig=P32081_PF00313_singleDEL$FDR<=.05
  P32081_PF00313_singleDEL<-P32081_PF00313_singleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #doubleDEL
  zstats=(P32081_PF00313_doubleDEL$scaled_fitness-1)/P32081_PF00313_doubleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P32081_PF00313_doubleDEL$zstat=zstats
  P32081_PF00313_doubleDEL$pvalue=pvals_man
  P32081_PF00313_doubleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P32081_PF00313_doubleDEL$fdrsig=P32081_PF00313_doubleDEL$FDR<=.05
  P32081_PF00313_doubleDEL<-P32081_PF00313_doubleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #tripleDEL
  zstats=(P32081_PF00313_tripleDEL$scaled_fitness-1)/P32081_PF00313_tripleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P32081_PF00313_tripleDEL$zstat=zstats
  P32081_PF00313_tripleDEL$pvalue=pvals_man
  P32081_PF00313_tripleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P32081_PF00313_tripleDEL$fdrsig=P32081_PF00313_tripleDEL$FDR<=.05
  P32081_PF00313_tripleDEL<-P32081_PF00313_tripleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #DelSub
  zstats=(P32081_PF00313_Delsub$scaled_fitness-1)/P32081_PF00313_Delsub$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P32081_PF00313_Delsub$zstat=zstats
  P32081_PF00313_Delsub$pvalue=pvals_man
  P32081_PF00313_Delsub[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P32081_PF00313_Delsub$fdrsig=P32081_PF00313_Delsub$FDR<=.05
  P32081_PF00313_Delsub<-P32081_PF00313_Delsub %>% mutate(significant=
                                                            ifelse(fdrsig == "TRUE", "*",NA))
  
  ##some of the DelSubs can cause 2 different substitutions at the same positions. We need to seperate them out into 2 dataframes. 
  P32081_PF00313_Delsub_2<-P32081_PF00313_Delsub[duplicated(P32081_PF00313_Delsub$Pos),]
  P32081_PF00313_Delsub<-P32081_PF00313_Delsub[!duplicated(P32081_PF00313_Delsub$Pos),]
  P32081_PF00313_Delsub_2$type<-"delSub_1"
  P32081_PF00313_Delsub$type<-"delSub_2"
  P32081_PF00313_Delsub_2$Pos<-as.numeric(P32081_PF00313_Delsub_2$Pos)
  P32081_PF00313_Delsub$Pos<-as.numeric(P32081_PF00313_Delsub$Pos)
  
  P32081_PF00313_Delsub<-P32081_PF00313_Delsub[P32081_PF00313_Delsub$STOP==FALSE,]
  
  
  ## final df
  P32081_PF00313_allvariants<-rbind(P32081_PF00313_syns[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
                                    P32081_PF00313_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P32081_PF00313_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P32081_PF00313_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P32081_PF00313_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P32081_PF00313_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")], 
                                    P32081_PF00313_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P32081_PF00313_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P32081_PF00313_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    fill=TRUE)
  
  return(list(P32081_PF00313_allvariants = P32081_PF00313_allvariants, 
              P32081_PF00313_syns = P32081_PF00313_syns[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
              P32081_PF00313_CX = P32081_PF00313_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P32081_PF00313_CXX = P32081_PF00313_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P32081_PF00313_CXXX = P32081_PF00313_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P32081_PF00313_singleDEL = P32081_PF00313_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P32081_PF00313_doubleDEL = P32081_PF00313_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P32081_PF00313_tripleDEL = P32081_PF00313_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P32081_PF00313_Delsub = P32081_PF00313_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P32081_PF00313_Delsub_2 = P32081_PF00313_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")]))
  
}

  
  
  
  
  
  
  
  
  
  
  