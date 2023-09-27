normalise_P01053_PF00280_data <- function() {
  
  ##isolate only the programmed lenghts and split into different dfs depending on indel type
  P01053_PF00280_variants$lenght=str_length(P01053_PF00280_variants$aa_seq)
  
  P01053_PF00280_singleDEL <- P01053_PF00280_variants[P01053_PF00280_variants$lenght == "62",]
  P01053_PF00280_doubleDEL <- P01053_PF00280_variants[P01053_PF00280_variants$lenght == "61",]
  P01053_PF00280_tripleDEL <- P01053_PF00280_variants[P01053_PF00280_variants$lenght == "60",]
  P01053_PF00280_CX <- P01053_PF00280_variants[P01053_PF00280_variants$lenght == "64",]
  P01053_PF00280_CXX <- P01053_PF00280_variants[P01053_PF00280_variants$lenght == "65",]
  P01053_PF00280_CXXX <- P01053_PF00280_variants[P01053_PF00280_variants$lenght == "66",]
  
  P01053_PF00280_singleDEL<-P01053_PF00280_singleDEL[P01053_PF00280_singleDEL$characters<3,]
  
  ##Delsubs
  P01053_PF00280_Delsub<-P01053_PF00280_variants[P01053_PF00280_variants$lenght == "62",]
  P01053_PF00280_Delsub<-P01053_PF00280_Delsub[P01053_PF00280_Delsub$characters>2,]
  
  ##extract position of mutation, add Pos  and length column
  P01053_PF00280_syns$lenght<-str_length(P01053_PF00280_syns$aa_seq)
  P01053_PF00280_syns$diff<-as.character(P01053_PF00280_syns$diff)
  P01053_PF00280_syns$diff<-str_first_number(P01053_PF00280_syns$diff)
  P01053_PF00280_syns$diff<-as.numeric(P01053_PF00280_syns$diff)
  
  ##extract positions as numeric into the Pos column for all df
  P01053_PF00280_singleDEL$Pos<-as.numeric(P01053_PF00280_singleDEL$diff)
  
  P01053_PF00280_doubleDEL$diff<-as.character(P01053_PF00280_doubleDEL$diff)
  P01053_PF00280_doubleDEL$Pos<-str_first_number(P01053_PF00280_doubleDEL$diff)
  P01053_PF00280_doubleDEL$Pos<-as.numeric(P01053_PF00280_doubleDEL$Pos)
  
  P01053_PF00280_tripleDEL$diff<-as.character(P01053_PF00280_tripleDEL$diff)
  P01053_PF00280_tripleDEL$Pos<-str_first_number(P01053_PF00280_tripleDEL$diff)
  P01053_PF00280_tripleDEL$Pos<-as.numeric(P01053_PF00280_tripleDEL$Pos)
  
  P01053_PF00280_CX$Pos<-as.numeric(P01053_PF00280_CX$diff)
  
  P01053_PF00280_CXX$diff<-as.character(P01053_PF00280_CXX$diff)
  P01053_PF00280_CXX$Pos<-str_first_number(P01053_PF00280_CXX$diff)
  P01053_PF00280_CXX$Pos<-as.numeric(P01053_PF00280_CXX$Pos)
  
  P01053_PF00280_CXXX$diff<-as.character(P01053_PF00280_CXXX$diff)
  P01053_PF00280_CXXX$Pos<-str_first_number(P01053_PF00280_CXXX$diff)
  P01053_PF00280_CXXX$Pos<-as.numeric(P01053_PF00280_CXXX$Pos)

  P01053_PF00280_Delsub$diff<-as.character(P01053_PF00280_Delsub$diff)
  P01053_PF00280_Delsub$Pos<-str_first_number(P01053_PF00280_Delsub$diff)
  
  ### find the aa position of the synonymous mutation
  P01053_PF00280_syns$Pos <- as.character(floor(P01053_PF00280_syns$diff / 3))
  
  
  ## add mutation type
  P01053_PF00280_singleDEL$type<-"singleDEL"
  P01053_PF00280_doubleDEL$type<-"doubleDEL"
  P01053_PF00280_tripleDEL$type<-"tripleDEL"
  P01053_PF00280_syns$type<-"synonymous"
  P01053_PF00280_CX$type<-"singleINS"
  P01053_PF00280_CXX$type<-"doubleINS"
  P01053_PF00280_CXXX$type<-"tripleINS"
  P01053_PF00280_Delsub$type<-"delSub"    
  
  ##fix mutation positions that are wrong manually. 
  # this is just an artifact of the function used to call the positions of difference and differing aa sequence
  P01053_PF00280_CXX[P01053_PF00280_CXX$diff == "c(13, 15)"]$Pos<-14
  P01053_PF00280_CXX[P01053_PF00280_CXX$diff == "c(16, 18)"]$Pos<-17
  P01053_PF00280_CXX[P01053_PF00280_CXX$diff == "c(28, 30)"]$Pos<-29
  
  P01053_PF00280_CXXX[P01053_PF00280_CXXX$diff == "c(13, 15, 16)",]$Pos<-14
  P01053_PF00280_CXXX[P01053_PF00280_CXXX$diff == "c(16, 18, 19)",]$Pos<-17
  P01053_PF00280_CXXX[P01053_PF00280_CXXX$diff == "c(28, 30, 31)",]$Pos<-29
  P01053_PF00280_CXXX[P01053_PF00280_CXXX$diff == "c(45, 46, 48)",]$Pos<-46
  
  P01053_PF00280_doubleDEL[P01053_PF00280_doubleDEL$diff == "c(13, 15)",]$Pos<-14
  P01053_PF00280_doubleDEL[P01053_PF00280_doubleDEL$diff == "c(16, 18)",]$Pos<-17
  P01053_PF00280_doubleDEL[P01053_PF00280_doubleDEL$diff == "c(28, 30)",]$Pos<-29
  
  P01053_PF00280_tripleDEL[P01053_PF00280_tripleDEL$diff == "c(13, 15, 16)",]$Pos<-14
  P01053_PF00280_tripleDEL[P01053_PF00280_tripleDEL$diff == "c(16, 18, 19)",]$Pos<-17
  P01053_PF00280_tripleDEL[P01053_PF00280_tripleDEL$diff == "c(28, 30, 31)",]$Pos<-29
  P01053_PF00280_tripleDEL[P01053_PF00280_tripleDEL$diff == "c(45, 46, 48)",]$Pos<-46
  
  ## add missing variants to the df
  # these are variants that have the same aa_seq although the mutation is technically in different positons. 
  #singleDEL: 13=14, 16=17, 28=29, 
  P01053_PF00280_singleDEL_missing<-P01053_PF00280_singleDEL[c(54, 53, 47),]
  P01053_PF00280_singleDEL_missing$Pos<-c(14,17,29)
  P01053_PF00280_singleDEL<-rbind(P01053_PF00280_singleDEL,
                                  P01053_PF00280_singleDEL_missing)
  #doubleDEL: 45=46
  P01053_PF00280_doubleDEL_missing<-P01053_PF00280_doubleDEL[c(40),]
  P01053_PF00280_doubleDEL_missing$Pos<-c(46)
  P01053_PF00280_doubleDEL<-rbind(P01053_PF00280_doubleDEL,
                                  P01053_PF00280_doubleDEL_missing)
  
  #tripleDEL: 3=4, 30=31, 35=36, 42=43, 51=52, 59=60
  P01053_PF00280_tripleDEL_missing<-P01053_PF00280_tripleDEL[c(53, 12, 15, 36, 31, 27),]
  P01053_PF00280_tripleDEL_missing$Pos<-c(4,31,36,43,52,60)
  P01053_PF00280_tripleDEL<-rbind(P01053_PF00280_tripleDEL,
                                  P01053_PF00280_tripleDEL_missing)
  #CX 13=14, 16=17, 28=29,
  P01053_PF00280_CX_missing<-P01053_PF00280_CX[c(7, 8, 14),]
  P01053_PF00280_CX_missing$Pos<-c(14,17,29)
  P01053_PF00280_CX<-rbind(P01053_PF00280_CX,
                           P01053_PF00280_CX_missing)
  #CXX  45=46
  P01053_PF00280_CXX_missing<-P01053_PF00280_CXX[c(22),]
  P01053_PF00280_CXX_missing$Pos<-c(46)
  P01053_PF00280_CXX<-rbind(P01053_PF00280_CXX,
                            P01053_PF00280_CXX_missing)
  #CXXX 3=4, 30=31, 35=36, 42=43, 51:52, 59=60
  P01053_PF00280_CXXX_missing<-P01053_PF00280_CXXX[c(3, 44, 41, 20,25, 29),]
  P01053_PF00280_CXXX_missing$Pos<-c(4,31,36,43,52,60)
  P01053_PF00280_CXXX<-rbind(P01053_PF00280_CXXX,
                             P01053_PF00280_CXXX_missing)
  
  ## data normalisation
  # normalise the data to mode of dead peak of the bimodal distribution of variants
  test<-rbind(P01053_PF00280_syns,
              P01053_PF00280_CX,
              P01053_PF00280_CXX,
              P01053_PF00280_CXXX,
              P01053_PF00280_singleDEL,
              P01053_PF00280_doubleDEL, 
              P01053_PF00280_tripleDEL,
              fill=TRUE)
  
  ##find mode of the lower peak in the bimodal distribution --> referred to as STOPs
  STOPs<-mlv(na.omit(test$fitness), method="naive")
  
  ## substract the mode from all variants. 
  P01053_PF00280_CX$norm_fitness<-P01053_PF00280_CX$fitness-STOPs
  P01053_PF00280_CXX$norm_fitness<-P01053_PF00280_CXX$fitness-STOPs
  P01053_PF00280_CXXX$norm_fitness<-P01053_PF00280_CXXX$fitness-STOPs
  P01053_PF00280_singleDEL$norm_fitness<-P01053_PF00280_singleDEL$fitness-STOPs
  P01053_PF00280_doubleDEL$norm_fitness<-P01053_PF00280_doubleDEL$fitness-STOPs
  P01053_PF00280_tripleDEL$norm_fitness<-P01053_PF00280_tripleDEL$fitness-STOPs
  P01053_PF00280_Delsub$norm_fitness<-P01053_PF00280_Delsub$fitness-STOPs
  P01053_PF00280_syns$norm_fitness<-P01053_PF00280_syns$fitness-STOPs
  
  ## find weighted mean of the synonymous variants
  syns<-sum(P01053_PF00280_syns$norm_fitness/(P01053_PF00280_syns$sigma^2))/sum(1/(P01053_PF00280_syns$sigma^2))
  
  ##normalise by dividing variants with the weighted mean of the synonymous variants
  P01053_PF00280_syns$scaled_fitness<-P01053_PF00280_syns$norm_fitness/sqrt(syns)
  P01053_PF00280_CX$scaled_fitness<-P01053_PF00280_CX$norm_fitness/sqrt(syns)
  P01053_PF00280_CXX$scaled_fitness<-P01053_PF00280_CXX$norm_fitness/sqrt(syns)
  P01053_PF00280_CXXX$scaled_fitness<-P01053_PF00280_CXXX$norm_fitness/sqrt(syns)
  P01053_PF00280_singleDEL$scaled_fitness<-P01053_PF00280_singleDEL$norm_fitness/sqrt(syns)
  P01053_PF00280_doubleDEL$scaled_fitness<-P01053_PF00280_doubleDEL$norm_fitness/sqrt(syns)
  P01053_PF00280_tripleDEL$scaled_fitness<-P01053_PF00280_tripleDEL$norm_fitness/sqrt(syns)
  P01053_PF00280_Delsub$scaled_fitness<-P01053_PF00280_Delsub$norm_fitness/sqrt(syns)
  
  ## scale the error by dividing by the same factor: dividing by the squared root of the synonmous fitness
  P01053_PF00280_syns$scaled_sigma<-P01053_PF00280_syns$sigma/syns
  P01053_PF00280_CX$scaled_sigma<-P01053_PF00280_CX$sigma/syns
  P01053_PF00280_CXX$scaled_sigma<-P01053_PF00280_CXX$sigma/syns
  P01053_PF00280_CXXX$scaled_sigma<-P01053_PF00280_CXXX$sigma/syns
  P01053_PF00280_singleDEL$scaled_sigma<-P01053_PF00280_singleDEL$sigma/syns
  P01053_PF00280_doubleDEL$scaled_sigma<-P01053_PF00280_doubleDEL$sigma/syns
  P01053_PF00280_tripleDEL$scaled_sigma<-P01053_PF00280_tripleDEL$sigma/syns
  P01053_PF00280_Delsub$scaled_sigma<-P01053_PF00280_Delsub$sigma/syns
  
 
  #### test if the change in fitness is signifcant
  ## we need to test for significance for values > and < than 1.  
  ## We calculate a z-stat and then do a 2-tailed test
  # we do mulitple testing correction using bonferroni
  
  # CX
  zstats=(P01053_PF00280_CX$scaled_fitness-1)/P01053_PF00280_CX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P01053_PF00280_CX$zstat=zstats
  P01053_PF00280_CX$pvalue=pvals_man
  P01053_PF00280_CX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P01053_PF00280_CX$fdrsig=P01053_PF00280_CX$FDR<=.05
  P01053_PF00280_CX<-P01053_PF00280_CX %>% mutate(significant=
                                                    ifelse(fdrsig == "TRUE", "*",NA))
  
  
  ##CXX
  zstats=(P01053_PF00280_CXX$scaled_fitness-1)/P01053_PF00280_CXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P01053_PF00280_CXX$zstat=zstats
  P01053_PF00280_CXX$pvalue=pvals_man
  P01053_PF00280_CXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P01053_PF00280_CXX$fdrsig=P01053_PF00280_CXX$FDR<=.05
  P01053_PF00280_CXX<-P01053_PF00280_CXX %>% mutate(significant=
                                                      ifelse(fdrsig == "TRUE", "*",NA))
  ##CXXX
  zstats=(P01053_PF00280_CXXX$scaled_fitness-1)/P01053_PF00280_CXXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P01053_PF00280_CXXX$zstat=zstats
  P01053_PF00280_CXXX$pvalue=pvals_man
  P01053_PF00280_CXXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P01053_PF00280_CXXX$fdrsig=P01053_PF00280_CXXX$FDR<=.05
  P01053_PF00280_CXXX<-P01053_PF00280_CXXX %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #singleDEL
  zstats=(P01053_PF00280_singleDEL$scaled_fitness-1)/P01053_PF00280_singleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P01053_PF00280_singleDEL$zstat=zstats
  P01053_PF00280_singleDEL$pvalue=pvals_man
  P01053_PF00280_singleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P01053_PF00280_singleDEL$fdrsig=P01053_PF00280_singleDEL$FDR<=.05
  P01053_PF00280_singleDEL<-P01053_PF00280_singleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #doubleDEL
  zstats=(P01053_PF00280_doubleDEL$scaled_fitness-1)/P01053_PF00280_doubleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P01053_PF00280_doubleDEL$zstat=zstats
  P01053_PF00280_doubleDEL$pvalue=pvals_man
  P01053_PF00280_doubleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P01053_PF00280_doubleDEL$fdrsig=P01053_PF00280_doubleDEL$FDR<=.05
  P01053_PF00280_doubleDEL<-P01053_PF00280_doubleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #tripleDEL
  zstats=(P01053_PF00280_tripleDEL$scaled_fitness-1)/P01053_PF00280_tripleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P01053_PF00280_tripleDEL$zstat=zstats
  P01053_PF00280_tripleDEL$pvalue=pvals_man
  P01053_PF00280_tripleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P01053_PF00280_tripleDEL$fdrsig=P01053_PF00280_tripleDEL$FDR<=.05
  P01053_PF00280_tripleDEL<-P01053_PF00280_tripleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #DelSub
  zstats=(P01053_PF00280_Delsub$scaled_fitness-1)/P01053_PF00280_Delsub$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P01053_PF00280_Delsub$zstat=zstats
  P01053_PF00280_Delsub$pvalue=pvals_man
  P01053_PF00280_Delsub[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P01053_PF00280_Delsub$fdrsig=P01053_PF00280_Delsub$FDR<=.05
  P01053_PF00280_Delsub<-P01053_PF00280_Delsub %>% mutate(significant=
                                                            ifelse(fdrsig == "TRUE", "*",NA))
  
  
  
  P01053_PF00280_Delsub<-P01053_PF00280_Delsub[P01053_PF00280_Delsub$STOP==FALSE,]
  
  ## final df
  P01053_PF00280_allvariants<-rbind(P01053_PF00280_syns[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
                                    P01053_PF00280_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P01053_PF00280_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P01053_PF00280_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P01053_PF00280_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P01053_PF00280_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")], 
                                    P01053_PF00280_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P01053_PF00280_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    fill=TRUE)
  
  return(list(P01053_PF00280_allvariants = P01053_PF00280_allvariants, 
              P01053_PF00280_syns = P01053_PF00280_syns[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
              P01053_PF00280_CX = P01053_PF00280_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P01053_PF00280_CXX = P01053_PF00280_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P01053_PF00280_CXXX = P01053_PF00280_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P01053_PF00280_singleDEL = P01053_PF00280_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P01053_PF00280_doubleDEL = P01053_PF00280_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P01053_PF00280_tripleDEL = P01053_PF00280_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P01053_PF00280_Delsub = P01053_PF00280_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")]))
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  