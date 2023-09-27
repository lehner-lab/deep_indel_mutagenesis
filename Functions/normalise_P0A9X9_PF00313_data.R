normalise_P0A9X9_PF00313_data <- function() {
  
  ##isolate only the programmed lenghts and split into different dfs depending on indel type
  P0A9X9_PF00313_variants$lenght=str_length(P0A9X9_PF00313_variants$aa_seq)
  
  P0A9X9_PF00313_singleDEL <- P0A9X9_PF00313_variants[P0A9X9_PF00313_variants$lenght == "64",]
  P0A9X9_PF00313_doubleDEL <- P0A9X9_PF00313_variants[P0A9X9_PF00313_variants$lenght == "63",]
  P0A9X9_PF00313_tripleDEL <- P0A9X9_PF00313_variants[P0A9X9_PF00313_variants$lenght == "62",]
  P0A9X9_PF00313_CX <- P0A9X9_PF00313_variants[P0A9X9_PF00313_variants$lenght == "66",]
  P0A9X9_PF00313_CXX <- P0A9X9_PF00313_variants[P0A9X9_PF00313_variants$lenght == "67",]
  P0A9X9_PF00313_CXXX <- P0A9X9_PF00313_variants[P0A9X9_PF00313_variants$lenght == "68",]
  
  P0A9X9_PF00313_singleDEL<-P0A9X9_PF00313_singleDEL[P0A9X9_PF00313_singleDEL$characters<3,]
  
  ##Delsubs
  P0A9X9_PF00313_Delsub<-P0A9X9_PF00313_variants[P0A9X9_PF00313_variants$lenght == "64",]
  P0A9X9_PF00313_Delsub<-P0A9X9_PF00313_Delsub[P0A9X9_PF00313_Delsub$characters>2,]
  
  ##extract position of mutation, add Pos  and length column
  P0A9X9_PF00313_syns$lenght<-str_length(P0A9X9_PF00313_syns$aa_seq)
  P0A9X9_PF00313_syns$diff<-as.character(P0A9X9_PF00313_syns$diff)
  P0A9X9_PF00313_syns$diff<-str_first_number(P0A9X9_PF00313_syns$diff)
  P0A9X9_PF00313_syns$diff<-as.numeric(P0A9X9_PF00313_syns$diff)
  P0A9X9_PF00313_syns$Pos<-as.numeric(P0A9X9_PF00313_syns$diff)
  
  ##extract positions as numeric into the Pos column for all df
  P0A9X9_PF00313_singleDEL$Pos<-as.numeric(P0A9X9_PF00313_singleDEL$diff)
  
  P0A9X9_PF00313_doubleDEL$diff<-as.character(P0A9X9_PF00313_doubleDEL$diff)
  P0A9X9_PF00313_doubleDEL$Pos<-str_first_number(P0A9X9_PF00313_doubleDEL$diff)
  P0A9X9_PF00313_doubleDEL$Pos<-as.numeric(P0A9X9_PF00313_doubleDEL$Pos)
  
  P0A9X9_PF00313_tripleDEL$diff<-as.character(P0A9X9_PF00313_tripleDEL$diff)
  P0A9X9_PF00313_tripleDEL$Pos<-str_first_number(P0A9X9_PF00313_tripleDEL$diff)
  P0A9X9_PF00313_tripleDEL$Pos<-as.numeric(P0A9X9_PF00313_tripleDEL$Pos)

  P0A9X9_PF00313_CX$Pos<-as.numeric(P0A9X9_PF00313_CX$diff)
  
  P0A9X9_PF00313_CXX$diff<-as.character(P0A9X9_PF00313_CXX$diff)
  P0A9X9_PF00313_CXX$Pos<-str_first_number(P0A9X9_PF00313_CXX$diff)
  P0A9X9_PF00313_CXX$Pos<-as.numeric(P0A9X9_PF00313_CXX$Pos)
  
  P0A9X9_PF00313_CXXX$diff<-as.character(P0A9X9_PF00313_CXXX$diff)
  P0A9X9_PF00313_CXXX$Pos<-str_first_number(P0A9X9_PF00313_CXXX$diff)
  P0A9X9_PF00313_CXXX$Pos<-as.numeric(P0A9X9_PF00313_CXXX$Pos)
  
  P0A9X9_PF00313_Delsub$diff<-as.character(P0A9X9_PF00313_Delsub$diff)
  P0A9X9_PF00313_Delsub$Pos<-str_first_number(P0A9X9_PF00313_Delsub$diff)
  
  ### find the aa position of the synonymous mutation
  P0A9X9_PF00313_syns$Pos <- as.character(floor(P0A9X9_PF00313_syns$diff / 3))
  
  
  ## add mutation type
  P0A9X9_PF00313_singleDEL$type<-"singleDEL"
  P0A9X9_PF00313_doubleDEL$type<-"doubleDEL"
  P0A9X9_PF00313_tripleDEL$type<-"tripleDEL"
  P0A9X9_PF00313_syns$type<-"synonymous"
  P0A9X9_PF00313_CX$type<-"singleINS"
  P0A9X9_PF00313_CXX$type<-"doubleINS"
  P0A9X9_PF00313_CXXX$type<-"tripleINS"
  P0A9X9_PF00313_Delsub$type<-"delSub"   
  
  ##fix mutation positions that are wrong manually. 
  # this is just an artifact of the function used to call the positions of difference and differing aa sequence
  P0A9X9_PF00313_CXX[P0A9X9_PF00313_CXX$diff == "c(20, 22)",]$Pos<-21
  P0A9X9_PF00313_CXX[P0A9X9_PF00313_CXX$diff == "c(59, 61)",]$Pos<-60
  
  P0A9X9_PF00313_CXXX[P0A9X9_PF00313_CXXX$diff == "c(13, 14, 17)",]$Pos<-14
  P0A9X9_PF00313_CXXX[P0A9X9_PF00313_CXXX$diff == "c(13, 14, 16)",]$Pos<-15
  P0A9X9_PF00313_CXXX[P0A9X9_PF00313_CXXX$diff == "c(20, 22, 23)",]$Pos<-21
  P0A9X9_PF00313_CXXX[P0A9X9_PF00313_CXXX$diff == "c(26, 27, 29)",]$Pos<-27
  P0A9X9_PF00313_CXXX[P0A9X9_PF00313_CXXX$diff == "c(59, 61, 62)",]$Pos<-60
  
  P0A9X9_PF00313_doubleDEL[P0A9X9_PF00313_doubleDEL$diff == "c(20, 22)",]$Pos<-21
  P0A9X9_PF00313_doubleDEL[P0A9X9_PF00313_doubleDEL$diff == "c(59, 61)",]$Pos<-60
  
  P0A9X9_PF00313_tripleDEL[P0A9X9_PF00313_tripleDEL$diff == "c(13, 14, 16)",]$Pos<-14
  P0A9X9_PF00313_tripleDEL[P0A9X9_PF00313_tripleDEL$diff == "c(13, 14, 17)",]$Pos<-15
  P0A9X9_PF00313_tripleDEL[P0A9X9_PF00313_tripleDEL$diff == "c(20, 22, 23)",]$Pos<-21
  P0A9X9_PF00313_tripleDEL[P0A9X9_PF00313_tripleDEL$diff == "c(26, 27, 29)",]$Pos<-27
  P0A9X9_PF00313_tripleDEL[P0A9X9_PF00313_tripleDEL$diff == "c(59, 61, 62)",]$Pos<-60
  
  ## add missing variants to the df
  # these are variants that have the same aa_seq although the mutation is technically in different positons. 
  #singleDEL: 20=21, 59=60
  P0A9X9_PF00313_singleDEL_missing<-P0A9X9_PF00313_singleDEL[c(53, 32),]
  P0A9X9_PF00313_singleDEL_missing$Pos<-c(21,60)
  P0A9X9_PF00313_singleDEL<-rbind(P0A9X9_PF00313_singleDEL,
                                  P0A9X9_PF00313_singleDEL_missing)
  #doubleDEL:13=14, 13=15, 26=27
  P0A9X9_PF00313_doubleDEL_missing<-P0A9X9_PF00313_doubleDEL[c(5, 5, 9),]
  P0A9X9_PF00313_doubleDEL_missing$Pos<-c(14,15,27)
  P0A9X9_PF00313_doubleDEL<-rbind(P0A9X9_PF00313_doubleDEL,
                                  P0A9X9_PF00313_doubleDEL_missing)
  #tripleDEL: 27=28, 54=55
  P0A9X9_PF00313_tripleDEL_missing<-P0A9X9_PF00313_tripleDEL[c(46,24),]
  P0A9X9_PF00313_tripleDEL_missing$Pos<-c(28,55)
  P0A9X9_PF00313_tripleDEL<-rbind(P0A9X9_PF00313_tripleDEL,
                                  P0A9X9_PF00313_tripleDEL_missing)
  #CX 20=21, 59=60
  P0A9X9_PF00313_CX_missing<-P0A9X9_PF00313_CX[c(11, 32),]
  P0A9X9_PF00313_CX_missing$Pos<-c(21,60)
  P0A9X9_PF00313_CX<-rbind(P0A9X9_PF00313_CX,
                           P0A9X9_PF00313_CX_missing)
  #CXX  13=14, 13=15, 26=27
  P0A9X9_PF00313_CXX_missing<-P0A9X9_PF00313_CXX[c(57, 57, 53),]
  P0A9X9_PF00313_CXX_missing$Pos<-c(14,15,27)
  P0A9X9_PF00313_CXX<-rbind(P0A9X9_PF00313_CXX,
                            P0A9X9_PF00313_CXX_missing)
  #CXXX 27=28, 54=55
  P0A9X9_PF00313_CXXX_missing<-P0A9X9_PF00313_CXXX[c(16,38),]
  P0A9X9_PF00313_CXXX_missing$Pos<-c(28,55)
  P0A9X9_PF00313_CXXX<-rbind(P0A9X9_PF00313_CXXX,
                             P0A9X9_PF00313_CXXX_missing)
  
  ## data normalisation
  # normalise the data to mode of dead peak of the bimodal distribution of variants
  test<-rbind(P0A9X9_PF00313_syns,
              P0A9X9_PF00313_CX,
              P0A9X9_PF00313_CXX,
              P0A9X9_PF00313_CXXX,
              P0A9X9_PF00313_singleDEL,
              P0A9X9_PF00313_doubleDEL, 
              P0A9X9_PF00313_tripleDEL,
              fill=TRUE)
  
  ##find mode of the lower peak in the bimodal distribution --> referred to as STOPs
  STOPs<-mlv(na.omit(test$fitness), method="naive")
  
  ## substract the mode from all variants. 
  P0A9X9_PF00313_CX$norm_fitness<-P0A9X9_PF00313_CX$fitness-STOPs
  P0A9X9_PF00313_CXX$norm_fitness<-P0A9X9_PF00313_CXX$fitness-STOPs
  P0A9X9_PF00313_CXXX$norm_fitness<-P0A9X9_PF00313_CXXX$fitness-STOPs
  P0A9X9_PF00313_singleDEL$norm_fitness<-P0A9X9_PF00313_singleDEL$fitness-STOPs
  P0A9X9_PF00313_doubleDEL$norm_fitness<-P0A9X9_PF00313_doubleDEL$fitness-STOPs
  P0A9X9_PF00313_tripleDEL$norm_fitness<-P0A9X9_PF00313_tripleDEL$fitness-STOPs
  P0A9X9_PF00313_Delsub$norm_fitness<-P0A9X9_PF00313_Delsub$fitness-STOPs
  P0A9X9_PF00313_syns$norm_fitness<-P0A9X9_PF00313_syns$fitness-STOPs
  
  ## find weighted mean of the synonymous variants
  syns<-sum(P0A9X9_PF00313_syns$norm_fitness/(P0A9X9_PF00313_syns$sigma^2))/sum(1/(P0A9X9_PF00313_syns$sigma^2))
  
  ##normalise by dividing variants with the weighted mean of the synonymous variants
  P0A9X9_PF00313_syns$scaled_fitness<-P0A9X9_PF00313_syns$norm_fitness/syns
  P0A9X9_PF00313_CX$scaled_fitness<-P0A9X9_PF00313_CX$norm_fitness/syns
  P0A9X9_PF00313_CXX$scaled_fitness<-P0A9X9_PF00313_CXX$norm_fitness/syns
  P0A9X9_PF00313_CXXX$scaled_fitness<-P0A9X9_PF00313_CXXX$norm_fitness/syns
  P0A9X9_PF00313_singleDEL$scaled_fitness<-P0A9X9_PF00313_singleDEL$norm_fitness/syns
  P0A9X9_PF00313_doubleDEL$scaled_fitness<-P0A9X9_PF00313_doubleDEL$norm_fitness/syns
  P0A9X9_PF00313_tripleDEL$scaled_fitness<-P0A9X9_PF00313_tripleDEL$norm_fitness/syns
  P0A9X9_PF00313_Delsub$scaled_fitness<-P0A9X9_PF00313_Delsub$norm_fitness/syns
  
  ## scale the error by dividing by the same factor: dividing by the squared root of the synonmous fitness
  P0A9X9_PF00313_syns$scaled_sigma<-P0A9X9_PF00313_syns$sigma/sqrt(syns)
  P0A9X9_PF00313_CX$scaled_sigma<-P0A9X9_PF00313_CX$sigma/sqrt(syns)
  P0A9X9_PF00313_CXX$scaled_sigma<-P0A9X9_PF00313_CXX$sigma/sqrt(syns)
  P0A9X9_PF00313_CXXX$scaled_sigma<-P0A9X9_PF00313_CXXX$sigma/sqrt(syns)
  P0A9X9_PF00313_singleDEL$scaled_sigma<-P0A9X9_PF00313_singleDEL$sigma/sqrt(syns)
  P0A9X9_PF00313_doubleDEL$scaled_sigma<-P0A9X9_PF00313_doubleDEL$sigma/sqrt(syns)
  P0A9X9_PF00313_tripleDEL$scaled_sigma<-P0A9X9_PF00313_tripleDEL$sigma/sqrt(syns)
  P0A9X9_PF00313_Delsub$scaled_sigma<-P0A9X9_PF00313_Delsub$sigma/sqrt(syns)
  
  #### test if the change in fitness is signifcant
  ## we need to test for significance for values > and < than 1.  
  ## We calculate a z-stat and then do a 2-tailed test
  # we do mulitple testing correction using bonferroni

  
  # CX
  zstats=(P0A9X9_PF00313_CX$scaled_fitness-1)/P0A9X9_PF00313_CX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P0A9X9_PF00313_CX$zstat=zstats
  P0A9X9_PF00313_CX$pvalue=pvals_man
  P0A9X9_PF00313_CX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P0A9X9_PF00313_CX$fdrsig=P0A9X9_PF00313_CX$FDR<=.05
  P0A9X9_PF00313_CX<-P0A9X9_PF00313_CX %>% mutate(significant=
                                                    ifelse(fdrsig == "TRUE", "*",NA))
  
  
  ##CXX
  zstats=(P0A9X9_PF00313_CXX$scaled_fitness-1)/P0A9X9_PF00313_CXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P0A9X9_PF00313_CXX$zstat=zstats
  P0A9X9_PF00313_CXX$pvalue=pvals_man
  P0A9X9_PF00313_CXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P0A9X9_PF00313_CXX$fdrsig=P0A9X9_PF00313_CXX$FDR<=.05
  P0A9X9_PF00313_CXX<-P0A9X9_PF00313_CXX %>% mutate(significant=
                                                      ifelse(fdrsig == "TRUE", "*",NA))
  ##CXXX
  zstats=(P0A9X9_PF00313_CXXX$scaled_fitness-1)/P0A9X9_PF00313_CXXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P0A9X9_PF00313_CXXX$zstat=zstats
  P0A9X9_PF00313_CXXX$pvalue=pvals_man
  P0A9X9_PF00313_CXXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P0A9X9_PF00313_CXXX$fdrsig=P0A9X9_PF00313_CXXX$FDR<=.05
  P0A9X9_PF00313_CXXX<-P0A9X9_PF00313_CXXX %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #singleDEL
  zstats=(P0A9X9_PF00313_singleDEL$scaled_fitness-1)/P0A9X9_PF00313_singleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P0A9X9_PF00313_singleDEL$zstat=zstats
  P0A9X9_PF00313_singleDEL$pvalue=pvals_man
  P0A9X9_PF00313_singleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P0A9X9_PF00313_singleDEL$fdrsig=P0A9X9_PF00313_singleDEL$FDR<=.05
  P0A9X9_PF00313_singleDEL<-P0A9X9_PF00313_singleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #doubleDEL
  zstats=(P0A9X9_PF00313_doubleDEL$scaled_fitness-1)/P0A9X9_PF00313_doubleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P0A9X9_PF00313_doubleDEL$zstat=zstats
  P0A9X9_PF00313_doubleDEL$pvalue=pvals_man
  P0A9X9_PF00313_doubleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P0A9X9_PF00313_doubleDEL$fdrsig=P0A9X9_PF00313_doubleDEL$FDR<=.05
  P0A9X9_PF00313_doubleDEL<-P0A9X9_PF00313_doubleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #tripleDEL
  zstats=(P0A9X9_PF00313_tripleDEL$scaled_fitness-1)/P0A9X9_PF00313_tripleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P0A9X9_PF00313_tripleDEL$zstat=zstats
  P0A9X9_PF00313_tripleDEL$pvalue=pvals_man
  P0A9X9_PF00313_tripleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P0A9X9_PF00313_tripleDEL$fdrsig=P0A9X9_PF00313_tripleDEL$FDR<=.05
  P0A9X9_PF00313_tripleDEL<-P0A9X9_PF00313_tripleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #DelSub
  zstats=(P0A9X9_PF00313_Delsub$scaled_fitness-1)/P0A9X9_PF00313_Delsub$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P0A9X9_PF00313_Delsub$zstat=zstats
  P0A9X9_PF00313_Delsub$pvalue=pvals_man
  P0A9X9_PF00313_Delsub[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P0A9X9_PF00313_Delsub$fdrsig=P0A9X9_PF00313_Delsub$FDR<=.05
  P0A9X9_PF00313_Delsub<-P0A9X9_PF00313_Delsub %>% mutate(significant=
                                                            ifelse(fdrsig == "TRUE", "*",NA))
  
  ##some of the DelSubs can cause 2 different substitutions at the same positions. We need to seperate them out into 2 dataframes. 
  P0A9X9_PF00313_Delsub_2<-P0A9X9_PF00313_Delsub[duplicated(P0A9X9_PF00313_Delsub$Pos),]
  P0A9X9_PF00313_Delsub<-P0A9X9_PF00313_Delsub[!duplicated(P0A9X9_PF00313_Delsub$Pos),]
  P0A9X9_PF00313_Delsub_2$type<-"delSub_1"
  P0A9X9_PF00313_Delsub$type<-"delSub_2"
  P0A9X9_PF00313_Delsub_2$Pos<-as.numeric(P0A9X9_PF00313_Delsub_2$Pos)
  P0A9X9_PF00313_Delsub$Pos<-as.numeric(P0A9X9_PF00313_Delsub$Pos)
  
  
  P0A9X9_PF00313_Delsub<-P0A9X9_PF00313_Delsub[P0A9X9_PF00313_Delsub$STOP==FALSE,]
  
  
  ## final df
  P0A9X9_PF00313_allvariants<-rbind(P0A9X9_PF00313_syns[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
                                    P0A9X9_PF00313_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P0A9X9_PF00313_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P0A9X9_PF00313_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P0A9X9_PF00313_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P0A9X9_PF00313_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")], 
                                    P0A9X9_PF00313_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P0A9X9_PF00313_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P0A9X9_PF00313_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    fill=TRUE)
  
  return(list(P0A9X9_PF00313_allvariants = P0A9X9_PF00313_allvariants, 
              P0A9X9_PF00313_syns = P0A9X9_PF00313_syns[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
              P0A9X9_PF00313_CX = P0A9X9_PF00313_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P0A9X9_PF00313_CXX = P0A9X9_PF00313_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P0A9X9_PF00313_CXXX = P0A9X9_PF00313_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P0A9X9_PF00313_singleDEL = P0A9X9_PF00313_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P0A9X9_PF00313_doubleDEL = P0A9X9_PF00313_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P0A9X9_PF00313_tripleDEL = P0A9X9_PF00313_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P0A9X9_PF00313_Delsub = P0A9X9_PF00313_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P0A9X9_PF00313_Delsub_2 = P0A9X9_PF00313_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")]))
  
}
  
  
  
  
  
  
  
  