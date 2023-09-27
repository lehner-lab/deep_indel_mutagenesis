normalise_P02417_PF01281_data <- function() {
  
  ##isolate only the programmed lenghts and split into different dfs depending on indel type
  P02417_PF01281_variants$lenght=str_length(P02417_PF01281_variants$aa_seq)
  
  P02417_PF01281_singleDEL <- P02417_PF01281_variants[P02417_PF01281_variants$lenght == "46",]
  P02417_PF01281_doubleDEL <- P02417_PF01281_variants[P02417_PF01281_variants$lenght == "45",]
  P02417_PF01281_tripleDEL <- P02417_PF01281_variants[P02417_PF01281_variants$lenght == "44",]
  P02417_PF01281_CX <- P02417_PF01281_variants[P02417_PF01281_variants$lenght == "48",]
  P02417_PF01281_CXX <- P02417_PF01281_variants[P02417_PF01281_variants$lenght == "49",]
  P02417_PF01281_CXXX <- P02417_PF01281_variants[P02417_PF01281_variants$lenght == "50",]
  
  P02417_PF01281_singleDEL<-P02417_PF01281_singleDEL[P02417_PF01281_singleDEL$characters<3,]
  
  ##Delsubs
  P02417_PF01281_Delsub<-P02417_PF01281_variants[P02417_PF01281_variants$lenght == "46",]
  P02417_PF01281_Delsub<-P02417_PF01281_Delsub[P02417_PF01281_Delsub$characters>2,]
  
  
  ##extract position of mutation, add Pos  and length column
  P02417_PF01281_syns$diff<-as.character(P02417_PF01281_syns$diff)
  P02417_PF01281_syns$diff<-str_first_number(P02417_PF01281_syns$diff)
  P02417_PF01281_syns$diff<-as.numeric(P02417_PF01281_syns$diff)
  P02417_PF01281_syns$Pos<-as.numeric(P02417_PF01281_syns$diff)
  P02417_PF01281_syns$lenght<-str_length(P02417_PF01281_syns$aa_seq)
  
  ##extract positions as numeric into the Pos column for all df
  P02417_PF01281_singleDEL$Pos<-as.numeric(P02417_PF01281_singleDEL$diff)
  
  P02417_PF01281_doubleDEL$diff<-as.character(P02417_PF01281_doubleDEL$diff)
  P02417_PF01281_doubleDEL$Pos<-str_first_number(P02417_PF01281_doubleDEL$diff)
  P02417_PF01281_doubleDEL$Pos<-as.numeric(P02417_PF01281_doubleDEL$Pos)
  
  P02417_PF01281_tripleDEL$diff<-as.character(P02417_PF01281_tripleDEL$diff)
  P02417_PF01281_tripleDEL$Pos<-str_first_number(P02417_PF01281_tripleDEL$diff)
  P02417_PF01281_tripleDEL$Pos<-as.numeric(P02417_PF01281_tripleDEL$Pos)
  
  P02417_PF01281_CX$Pos<-as.numeric(P02417_PF01281_CX$diff)
  
  P02417_PF01281_CXX$diff<-as.character(P02417_PF01281_CXX$diff)
  P02417_PF01281_CXX$Pos<-str_first_number(P02417_PF01281_CXX$diff)
  P02417_PF01281_CXX$Pos<-as.numeric(P02417_PF01281_CXX$Pos)
  
  P02417_PF01281_CXXX$diff<-as.character(P02417_PF01281_CXXX$diff)
  P02417_PF01281_CXXX$Pos<-str_first_number(P02417_PF01281_CXXX$diff)
  P02417_PF01281_CXXX$Pos<-as.numeric(P02417_PF01281_CXXX$Pos)
  
  P02417_PF01281_Delsub$diff<-as.character(P02417_PF01281_Delsub$diff)
  P02417_PF01281_Delsub$Pos<-str_first_number(P02417_PF01281_Delsub$diff)
  
  ### find the aa position of the synonymous mutation
  P02417_PF01281_syns$Pos <- as.character(floor(P02417_PF01281_syns$diff / 3))
  
  
  ## add mutation type
  P02417_PF01281_singleDEL$type<-"singleDEL"
  P02417_PF01281_doubleDEL$type<-"doubleDEL"
  P02417_PF01281_tripleDEL$type<-"tripleDEL"
  P02417_PF01281_CX$type<-"singleINS"
  P02417_PF01281_CXX$type<-"doubleINS"
  P02417_PF01281_CXXX$type<-"tripleINS"
  P02417_PF01281_Delsub$type<-"delSub"    
  P02417_PF01281_syns$type<-"synonymous"
  
  ##fix mutation positions that are called wrong manually. 
  # this is just an artifact of the function used to call the positions of difference and differing aa sequence
  P02417_PF01281_CXX[P02417_PF01281_CXX$diff == "c(14, 16)",]$Pos<-15
  P02417_PF01281_CXX[P02417_PF01281_CXX$diff == "c(27, 29)",]$Pos<-28
  
  P02417_PF01281_CXXX[P02417_PF01281_CXXX$diff == "c(10, 11, 13)",]$Pos<-11
  P02417_PF01281_CXXX[P02417_PF01281_CXXX$diff == "c(10, 11, 14)",]$Pos<-12
  P02417_PF01281_CXXX[P02417_PF01281_CXXX$diff == "c(14, 16, 17)",]$Pos<-15
  P02417_PF01281_CXXX[P02417_PF01281_CXXX$diff == "c(27, 29, 30)",]$Pos<-28
  P02417_PF01281_CXXX[P02417_PF01281_CXXX$diff == "c(29, 30, 32)",]$Pos<-30
  
  P02417_PF01281_doubleDEL[P02417_PF01281_doubleDEL$diff == "c(14, 16)",]$Pos<-15
  P02417_PF01281_doubleDEL[P02417_PF01281_doubleDEL$diff == "c(27, 29)",]$Pos<-28
  
  P02417_PF01281_tripleDEL[P02417_PF01281_tripleDEL$diff == "c(10, 11, 13)",]$Pos<-11
  P02417_PF01281_tripleDEL[P02417_PF01281_tripleDEL$diff == "c(10, 11, 14)",]$Pos<-12
  P02417_PF01281_tripleDEL[P02417_PF01281_tripleDEL$diff == "c(14, 16, 17)",]$Pos<-15
  P02417_PF01281_tripleDEL[P02417_PF01281_tripleDEL$diff == "c(27, 29, 30)",]$Pos<-28
  P02417_PF01281_tripleDEL[P02417_PF01281_tripleDEL$diff == "c(29, 30, 32)",]$Pos<-30
  
  ## add missing variants to the df
  # these are variants that have the same aa_seq although the mutation is technically in different positons. 
  #singleDEL: 14=15, 27=28
  P02417_PF01281_singleDEL_missing<-P02417_PF01281_singleDEL[c(39, 34),]
  P02417_PF01281_singleDEL_missing$Pos<-c(15,28)
  P02417_PF01281_singleDEL<-rbind(P02417_PF01281_singleDEL,
                                  P02417_PF01281_singleDEL_missing)
  #doubleDEL:10=11, 10=12, 10=13, 29=30
  P02417_PF01281_doubleDEL_missing<-P02417_PF01281_doubleDEL[c(5,5,5,13),]
  P02417_PF01281_doubleDEL_missing$Pos<-c(11,12,13,30)
  P02417_PF01281_doubleDEL<-rbind(P02417_PF01281_doubleDEL,
                                  P02417_PF01281_doubleDEL_missing)
  #tripleDEL: 7=8, 12=13, 12=14, 36=37, 39=40, 44=45
  P02417_PF01281_tripleDEL_missing<-P02417_PF01281_tripleDEL[c(36, 34,34, 16, 18, 20),]
  P02417_PF01281_tripleDEL_missing$Pos<-c(8,13,14,37,40,45)
  P02417_PF01281_tripleDEL<-rbind(P02417_PF01281_tripleDEL,
                                  P02417_PF01281_tripleDEL_missing)
  #CX 14=15, 27=28 
  P02417_PF01281_CX_missing<-P02417_PF01281_CX[c(7, 12),]
  P02417_PF01281_CX_missing$Pos<-c(15,28)
  P02417_PF01281_CX<-rbind(P02417_PF01281_CX,
                           P02417_PF01281_CX_missing)
  #CXX 10=11, 10=12, 10=13, 29=30
  P02417_PF01281_CXX_missing<-P02417_PF01281_CXX[c(38,38,38,30),]
  P02417_PF01281_CXX_missing$Pos<-c(11,12,13,30)
  P02417_PF01281_CXX<-rbind(P02417_PF01281_CXX,
                            P02417_PF01281_CXX_missing)
  #CXXX 7=8, 12=13, 12=14, 36=37, 39=40, 44=45
  P02417_PF01281_CXXX_missing<-P02417_PF01281_CXXX[c(4,6,6,24, 22, 20),]
  P02417_PF01281_CXXX_missing$Pos<-c(8,13,14,37,40,45)
  P02417_PF01281_CXXX<-rbind(P02417_PF01281_CXXX,
                             P02417_PF01281_CXXX_missing)
  
  ## data normalisation
  # normalise the data to mode of dead peak of the bimodal distribution of variants
  test<-rbind(P02417_PF01281_syns,
              P02417_PF01281_CX,
              P02417_PF01281_CXX,
              P02417_PF01281_CXXX,
              P02417_PF01281_singleDEL,
              P02417_PF01281_doubleDEL, 
              P02417_PF01281_tripleDEL,
              fill=TRUE)
  
  
  #bind all data
  STOPs<-min(mlv(na.omit(test$fitness), method="naive"))
  
  ## substract the mode from all variants. 
  P02417_PF01281_CX$norm_fitness<-P02417_PF01281_CX$fitness-STOPs
  P02417_PF01281_CXX$norm_fitness<-P02417_PF01281_CXX$fitness-STOPs
  P02417_PF01281_CXXX$norm_fitness<-P02417_PF01281_CXXX$fitness-STOPs
  P02417_PF01281_singleDEL$norm_fitness<-P02417_PF01281_singleDEL$fitness-STOPs
  P02417_PF01281_doubleDEL$norm_fitness<-P02417_PF01281_doubleDEL$fitness-STOPs
  P02417_PF01281_tripleDEL$norm_fitness<-P02417_PF01281_tripleDEL$fitness-STOPs
  P02417_PF01281_Delsub$norm_fitness<-P02417_PF01281_Delsub$fitness-STOPs
  P02417_PF01281_syns$norm_fitness<-P02417_PF01281_syns$fitness-STOPs
  
  ## find weighted mean of the synonymous variants
  syns<-sum(P02417_PF01281_syns$norm_fitness/(P02417_PF01281_syns$sigma^2))/sum(1/(P02417_PF01281_syns$sigma^2))
  
  ##normalise by dividing variants with the weighted mean of the synonymous variants
  P02417_PF01281_syns$scaled_fitness<-P02417_PF01281_syns$norm_fitness/syns
  P02417_PF01281_CX$scaled_fitness<-P02417_PF01281_CX$norm_fitness/syns
  P02417_PF01281_CXX$scaled_fitness<-P02417_PF01281_CXX$norm_fitness/syns
  P02417_PF01281_CXXX$scaled_fitness<-P02417_PF01281_CXXX$norm_fitness/syns
  P02417_PF01281_singleDEL$scaled_fitness<-P02417_PF01281_singleDEL$norm_fitness/syns
  P02417_PF01281_doubleDEL$scaled_fitness<-P02417_PF01281_doubleDEL$norm_fitness/syns
  P02417_PF01281_tripleDEL$scaled_fitness<-P02417_PF01281_tripleDEL$norm_fitness/syns
  P02417_PF01281_Delsub$scaled_fitness<-P02417_PF01281_Delsub$norm_fitness/syns
  
  
  ## scale the error by dividing by the same factor: dividing by the squared root of the synonmous fitness
  P02417_PF01281_syns$scaled_sigma<-P02417_PF01281_syns$sigma/sqrt(syns)
  P02417_PF01281_CX$scaled_sigma<-P02417_PF01281_CX$sigma/sqrt(syns)
  P02417_PF01281_CXX$scaled_sigma<-P02417_PF01281_CXX$sigma/sqrt(syns)
  P02417_PF01281_CXXX$scaled_sigma<-P02417_PF01281_CXXX$sigma/sqrt(syns)
  P02417_PF01281_singleDEL$scaled_sigma<-P02417_PF01281_singleDEL$sigma/sqrt(syns)
  P02417_PF01281_doubleDEL$scaled_sigma<-P02417_PF01281_doubleDEL$sigma/sqrt(syns)
  P02417_PF01281_tripleDEL$scaled_sigma<-P02417_PF01281_tripleDEL$sigma/sqrt(syns)
  P02417_PF01281_Delsub$scaled_sigma<-P02417_PF01281_Delsub$sigma/sqrt(syns)
  
  #### test if the change in fitness is signifcant
  ## we need to test for significance for values > and < than 1.  
  ## We calculate a z-stat and then do a 2-tailed test
  # we do mulitple testing correction using bonferroni
  
  # CX
  zstats=(P02417_PF01281_CX$scaled_fitness-1)/P02417_PF01281_CX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P02417_PF01281_CX$zstat=zstats
  P02417_PF01281_CX$pvalue=pvals_man
  P02417_PF01281_CX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P02417_PF01281_CX$fdrsig=P02417_PF01281_CX$FDR<=.05
  P02417_PF01281_CX<-P02417_PF01281_CX %>% mutate(significant=
                                                    ifelse(fdrsig == "TRUE", "*",NA))
  
  
  ##CXX
  zstats=(P02417_PF01281_CXX$scaled_fitness-1)/P02417_PF01281_CXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P02417_PF01281_CXX$zstat=zstats
  P02417_PF01281_CXX$pvalue=pvals_man
  P02417_PF01281_CXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P02417_PF01281_CXX$fdrsig=P02417_PF01281_CXX$FDR<=.05
  P02417_PF01281_CXX<-P02417_PF01281_CXX %>% mutate(significant=
                                                      ifelse(fdrsig == "TRUE", "*",NA))
  ##CXXX
  zstats=(P02417_PF01281_CXXX$scaled_fitness-1)/P02417_PF01281_CXXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P02417_PF01281_CXXX$zstat=zstats
  P02417_PF01281_CXXX$pvalue=pvals_man
  P02417_PF01281_CXXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P02417_PF01281_CXXX$fdrsig=P02417_PF01281_CXXX$FDR<=.05
  P02417_PF01281_CXXX<-P02417_PF01281_CXXX %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #singleDEL
  zstats=(P02417_PF01281_singleDEL$scaled_fitness-1)/P02417_PF01281_singleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P02417_PF01281_singleDEL$zstat=zstats
  P02417_PF01281_singleDEL$pvalue=pvals_man
  P02417_PF01281_singleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P02417_PF01281_singleDEL$fdrsig=P02417_PF01281_singleDEL$FDR<=.05
  P02417_PF01281_singleDEL<-P02417_PF01281_singleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #doubleDEL
  zstats=(P02417_PF01281_doubleDEL$scaled_fitness-1)/P02417_PF01281_doubleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P02417_PF01281_doubleDEL$zstat=zstats
  P02417_PF01281_doubleDEL$pvalue=pvals_man
  P02417_PF01281_doubleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P02417_PF01281_doubleDEL$fdrsig=P02417_PF01281_doubleDEL$FDR<=.05
  P02417_PF01281_doubleDEL<-P02417_PF01281_doubleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #tripleDEL
  zstats=(P02417_PF01281_tripleDEL$scaled_fitness-1)/P02417_PF01281_tripleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P02417_PF01281_tripleDEL$zstat=zstats
  P02417_PF01281_tripleDEL$pvalue=pvals_man
  P02417_PF01281_tripleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P02417_PF01281_tripleDEL$fdrsig=P02417_PF01281_tripleDEL$FDR<=.05
  P02417_PF01281_tripleDEL<-P02417_PF01281_tripleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #DelSub
  zstats=(P02417_PF01281_Delsub$scaled_fitness-1)/P02417_PF01281_Delsub$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P02417_PF01281_Delsub$zstat=zstats
  P02417_PF01281_Delsub$pvalue=pvals_man
  P02417_PF01281_Delsub[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P02417_PF01281_Delsub$fdrsig=P02417_PF01281_Delsub$FDR<=.05
  P02417_PF01281_Delsub<-P02417_PF01281_Delsub %>% mutate(significant=
                                                            ifelse(fdrsig == "TRUE", "*",NA))
  
  ##some of the DelSubs can cause 2 different substitutions at the same positions. We need to seperate them out into 2 dataframes. 
  P02417_PF01281_Delsub_2<-P02417_PF01281_Delsub[duplicated(P02417_PF01281_Delsub$Pos),]
  P02417_PF01281_Delsub<-P02417_PF01281_Delsub[!duplicated(P02417_PF01281_Delsub$Pos),]
  P02417_PF01281_Delsub_2$type<-"delSub_1"
  P02417_PF01281_Delsub$type<-"delSub_2"
  P02417_PF01281_Delsub_2$Pos<-as.numeric(P02417_PF01281_Delsub_2$Pos)
  P02417_PF01281_Delsub$Pos<-as.numeric(P02417_PF01281_Delsub$Pos)
  
  P02417_PF01281_Delsub<-P02417_PF01281_Delsub[P02417_PF01281_Delsub$STOP==FALSE,]
  
  
  
  ## final df
  P02417_PF01281_allvariants<-rbind(P02417_PF01281_syns[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
                                    P02417_PF01281_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P02417_PF01281_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P02417_PF01281_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P02417_PF01281_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P02417_PF01281_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")], 
                                    P02417_PF01281_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P02417_PF01281_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P02417_PF01281_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    fill=TRUE)
  
  return(list(P02417_PF01281_allvariants = P02417_PF01281_allvariants, 
              P02417_PF01281_syns = P02417_PF01281_syns[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
              P02417_PF01281_CX = P02417_PF01281_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P02417_PF01281_CXX = P02417_PF01281_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P02417_PF01281_CXXX = P02417_PF01281_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P02417_PF01281_singleDEL = P02417_PF01281_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P02417_PF01281_doubleDEL = P02417_PF01281_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P02417_PF01281_tripleDEL = P02417_PF01281_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P02417_PF01281_Delsub = P02417_PF01281_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P02417_PF01281_Delsub_2 = P02417_PF01281_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")]))
  
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  