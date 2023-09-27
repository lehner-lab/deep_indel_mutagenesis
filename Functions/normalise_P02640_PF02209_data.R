normalise_P02640_PF02209_data <- function() {
  
  ##isolate only the programmed lenghts and split into different dfs depending on indel type
  P02640_PF02209_variants$lenght=str_length(P02640_PF02209_variants$aa_seq)
  
  P02640_PF02209_singleDEL <- P02640_PF02209_variants[P02640_PF02209_variants$lenght == "35",]
  P02640_PF02209_doubleDEL <- P02640_PF02209_variants[P02640_PF02209_variants$lenght == "34",]
  P02640_PF02209_tripleDEL <- P02640_PF02209_variants[P02640_PF02209_variants$lenght == "33",]
  P02640_PF02209_CX <- P02640_PF02209_variants[P02640_PF02209_variants$lenght == "37",]
  P02640_PF02209_CXX <- P02640_PF02209_variants[P02640_PF02209_variants$lenght == "38",]
  P02640_PF02209_CXXX <- P02640_PF02209_variants[P02640_PF02209_variants$lenght == "39",]
  
  
  P02640_PF02209_singleDEL<-P02640_PF02209_singleDEL[P02640_PF02209_singleDEL$characters<3,]
  
  ##Delsubs
  P02640_PF02209_Delsub<-P02640_PF02209_variants[P02640_PF02209_variants$lenght == "35",]
  P02640_PF02209_Delsub<-P02640_PF02209_Delsub[P02640_PF02209_Delsub$characters>2,]
  
  ##extract position of mutation, add Pos  and length column
  P02640_PF02209_syns$diff<-as.character(P02640_PF02209_syns$diff)
  P02640_PF02209_syns$diff<-str_first_number(P02640_PF02209_syns$diff)
  P02640_PF02209_syns$diff<-as.numeric(P02640_PF02209_syns$diff)
  P02640_PF02209_syns$Pos<-as.numeric(P02640_PF02209_syns$diff)
  P02640_PF02209_syns$lenght<-str_length(P02640_PF02209_syns$aa_seq)
  
  ##extract positions as numeric into the Pos column for all df
  P02640_PF02209_singleDEL$Pos<-as.numeric(P02640_PF02209_singleDEL$diff)
  
  P02640_PF02209_doubleDEL$diff<-as.character(P02640_PF02209_doubleDEL$diff)
  P02640_PF02209_doubleDEL$Pos<-str_first_number(P02640_PF02209_doubleDEL$diff)
  P02640_PF02209_doubleDEL$Pos<-as.numeric(P02640_PF02209_doubleDEL$Pos)
  
  P02640_PF02209_tripleDEL$diff<-as.character(P02640_PF02209_tripleDEL$diff)
  P02640_PF02209_tripleDEL$Pos<-str_first_number(P02640_PF02209_tripleDEL$diff)
  P02640_PF02209_tripleDEL$Pos<-as.numeric(P02640_PF02209_tripleDEL$Pos)
  
  P02640_PF02209_CX$Pos<-as.numeric(P02640_PF02209_CX$diff)
  
  P02640_PF02209_CXX$diff<-as.character(P02640_PF02209_CXX$diff)
  P02640_PF02209_CXX$Pos<-str_first_number(P02640_PF02209_CXX$diff)
  P02640_PF02209_CXX$Pos<-as.numeric(P02640_PF02209_CXX$Pos)
  
  P02640_PF02209_CXXX$diff<-as.character(P02640_PF02209_CXXX$diff)
  P02640_PF02209_CXXX$Pos<-str_first_number(P02640_PF02209_CXXX$diff)
  P02640_PF02209_CXXX$Pos<-as.numeric(P02640_PF02209_CXXX$Pos)
  
  P02640_PF02209_Delsub$diff<-as.character(P02640_PF02209_Delsub$diff)
  P02640_PF02209_Delsub$Pos<-str_first_number(P02640_PF02209_Delsub$diff)
  
  ### find the aa position of the synonymous mutation
  P02640_PF02209_syns$Pos <- as.character(floor(P02640_PF02209_syns$diff / 3))
  
  
  ## add mutation type
  P02640_PF02209_singleDEL$type<-"singleDEL"
  P02640_PF02209_doubleDEL$type<-"doubleDEL"
  P02640_PF02209_tripleDEL$type<-"tripleDEL"
  P02640_PF02209_CX$type<-"singleINS"
  P02640_PF02209_CXX$type<-"doubleINS"
  P02640_PF02209_CXXX$type<-"tripleINS"
  P02640_PF02209_Delsub$type<-"delSub"    
  P02640_PF02209_syns$type<-"synonymous"
  
  ##fix mutation positions that are called wrong manually. 
  # this is just an artifact of the function used to call the positions of difference and differing aa sequence
  P02640_PF02209_CXX[P02640_PF02209_CXX$diff == "c(26, 28)",]$Pos<-27
  P02640_PF02209_CXX[P02640_PF02209_CXX$diff == "c(30, 32)",]$Pos<-31
  
  P02640_PF02209_CXXX[P02640_PF02209_CXXX$diff == "c(4, 5, 7)",]$Pos<-5
  P02640_PF02209_CXXX[P02640_PF02209_CXXX$diff == "c(17, 18, 20)",]$Pos<-18
  P02640_PF02209_CXXX[P02640_PF02209_CXXX$diff == "c(21, 22, 24)",]$Pos<-22
  P02640_PF02209_CXXX[P02640_PF02209_CXXX$diff == "c(26, 28, 29)",]$Pos<-27
  P02640_PF02209_CXXX[P02640_PF02209_CXXX$diff == "c(30, 32, 34)",]$Pos<-31
  
  P02640_PF02209_doubleDEL[P02640_PF02209_doubleDEL$diff == "c(26, 28)",]$Pos<-27
  P02640_PF02209_doubleDEL[P02640_PF02209_doubleDEL$diff == "c(30, 32)",]$Pos<-31
  
  P02640_PF02209_tripleDEL[P02640_PF02209_tripleDEL$diff == "c(4, 5, 7)",]$Pos<-5
  P02640_PF02209_tripleDEL[P02640_PF02209_tripleDEL$diff == "c(17, 18, 20)",]$Pos<-18
  P02640_PF02209_tripleDEL[P02640_PF02209_tripleDEL$diff == "c(21, 22, 24)",]$Pos<-22
  P02640_PF02209_tripleDEL[P02640_PF02209_tripleDEL$diff == "c(26, 28, 29)",]$Pos<-27
  P02640_PF02209_tripleDEL[P02640_PF02209_tripleDEL$diff == "c(30, 32, 34)",]$Pos<-31
  
  ## add missing variants to the df
  # these are variants that have the same aa_seq although the mutation is technically in different positons. 
  #singleDEL: 26=27, 30=31
  P02640_PF02209_singleDEL_missing<-P02640_PF02209_singleDEL[c(14, 21),]
  P02640_PF02209_singleDEL_missing$Pos<-c(27,31)
  P02640_PF02209_singleDEL<-rbind(P02640_PF02209_singleDEL,
                                  P02640_PF02209_singleDEL_missing)
  #doubleDEL: 4=5, 17=18, 21=22, 31=32
  P02640_PF02209_doubleDEL_missing<-P02640_PF02209_doubleDEL[c(30, 9, 22, 18),]
  P02640_PF02209_doubleDEL_missing$Pos<-c(5,18,22,32)
  P02640_PF02209_doubleDEL<-rbind(P02640_PF02209_doubleDEL,
                                  P02640_PF02209_doubleDEL_missing)
  #tripleDEL: 30=31
  P02640_PF02209_tripleDEL_missing<-P02640_PF02209_tripleDEL[c(22),]
  P02640_PF02209_tripleDEL_missing$Pos<-c(31)
  P02640_PF02209_tripleDEL<-rbind(P02640_PF02209_tripleDEL,
                                  P02640_PF02209_tripleDEL_missing)
  #CX 26=27, 30=31
  P02640_PF02209_CX_missing<-P02640_PF02209_CX[c(21, 14),]
  P02640_PF02209_CX_missing$Pos<-c(27,31)
  P02640_PF02209_CX<-rbind(P02640_PF02209_CX,
                           P02640_PF02209_CX_missing)
  #CXX 4=5, 17=18, 21=22, 31=32
  P02640_PF02209_CXX_missing<-P02640_PF02209_CXX[c(2, 23,10, 14),]
  P02640_PF02209_CXX_missing$Pos<-c(5,18,22,32)
  P02640_PF02209_CXX<-rbind(P02640_PF02209_CXX,
                            P02640_PF02209_CXX_missing)
  #CXXX 30=31
  P02640_PF02209_CXXX_missing<-P02640_PF02209_CXXX[c(12),]
  P02640_PF02209_CXXX_missing$Pos<-c(31)
  P02640_PF02209_CXXX<-rbind(P02640_PF02209_CXXX,
                             P02640_PF02209_CXXX_missing)
  
  P02640_PF02209_tripleDEL[P02640_PF02209_tripleDEL$aa_seq == "HLSDEDFKAVFGMTRSAFANLPLWKQQNLKKLF", ]$Pos<-"32"
  P02640_PF02209_CXXX[P02640_PF02209_CXXX$aa_seq == "HLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGEKGLF", ]$Pos<-"32"
  
  ## data normalisation
  # normalise the data to mode of dead peak of the bimodal distribution of variants
  #bind all data
  test<-rbind(P02640_PF02209_syns,
              P02640_PF02209_CX,
              P02640_PF02209_CXX,
              P02640_PF02209_CXXX,
              P02640_PF02209_singleDEL,
              P02640_PF02209_doubleDEL, 
              P02640_PF02209_tripleDEL,
              fill=TRUE)
  
  ##find mode of the lower peak in the bimodal distribution --> referred to as STOPs
  #STOPs<-mlv(na.omit(test$fitness), method="naive")
  
  ### the biomodality of the data is an issue when determining the mode of the deleterious peak
  # fist do a density function, find the median and exclude observation above the median to exclude the "neutral effect" peak
  density_values <- density(test$fitness)
  # now find the mode using mlv function and the "Vieu" method
  # the "Vieu" method finds the point where the estimated density function has a local maximum
  # it is defined as the value at which the kernel density derivative estimate is null
  STOPs<-mlv(na.omit(test[test$fitness<mean(density_values$x),]$fitness), method="Vieu")
  
  ## substract the mode from all variants. 
  P02640_PF02209_CX$norm_fitness<-P02640_PF02209_CX$fitness-STOPs
  P02640_PF02209_CXX$norm_fitness<-P02640_PF02209_CXX$fitness-STOPs
  P02640_PF02209_CXXX$norm_fitness<-P02640_PF02209_CXXX$fitness-STOPs
  P02640_PF02209_singleDEL$norm_fitness<-P02640_PF02209_singleDEL$fitness-STOPs
  P02640_PF02209_doubleDEL$norm_fitness<-P02640_PF02209_doubleDEL$fitness-STOPs
  P02640_PF02209_tripleDEL$norm_fitness<-P02640_PF02209_tripleDEL$fitness-STOPs
  P02640_PF02209_Delsub$norm_fitness<-P02640_PF02209_Delsub$fitness-STOPs
  P02640_PF02209_syns$norm_fitness<-P02640_PF02209_syns$fitness-STOPs
  
  ## find weighted mean of the synonymous variants
  syns<-sum(P02640_PF02209_syns$norm_fitness/(P02640_PF02209_syns$sigma^2))/sum(1/(P02640_PF02209_syns$sigma^2))
  
  ##normalise by dividing variants with the weighted mean of the synonymous variants
  P02640_PF02209_syns$scaled_fitness<-P02640_PF02209_syns$norm_fitness/syns
  P02640_PF02209_CX$scaled_fitness<-P02640_PF02209_CX$norm_fitness/syns
  P02640_PF02209_CXX$scaled_fitness<-P02640_PF02209_CXX$norm_fitness/syns
  P02640_PF02209_CXXX$scaled_fitness<-P02640_PF02209_CXXX$norm_fitness/syns
  P02640_PF02209_singleDEL$scaled_fitness<-P02640_PF02209_singleDEL$norm_fitness/syns
  P02640_PF02209_doubleDEL$scaled_fitness<-P02640_PF02209_doubleDEL$norm_fitness/syns
  P02640_PF02209_tripleDEL$scaled_fitness<-P02640_PF02209_tripleDEL$norm_fitness/syns
  P02640_PF02209_Delsub$scaled_fitness<-P02640_PF02209_Delsub$norm_fitness/syns
  
  
  ## scale the error by dividing by the same factor: dividing by the squared root of the synonmous fitness
  P02640_PF02209_syns$scaled_sigma<-P02640_PF02209_syns$sigma/sqrt(syns)
  P02640_PF02209_CX$scaled_sigma<-P02640_PF02209_CX$sigma/sqrt(syns)
  P02640_PF02209_CXX$scaled_sigma<-P02640_PF02209_CXX$sigma/sqrt(syns)
  P02640_PF02209_CXXX$scaled_sigma<-P02640_PF02209_CXXX$sigma/sqrt(syns)
  P02640_PF02209_singleDEL$scaled_sigma<-P02640_PF02209_singleDEL$sigma/sqrt(syns)
  P02640_PF02209_doubleDEL$scaled_sigma<-P02640_PF02209_doubleDEL$sigma/sqrt(syns)
  P02640_PF02209_tripleDEL$scaled_sigma<-P02640_PF02209_tripleDEL$sigma/sqrt(syns)
  P02640_PF02209_Delsub$scaled_sigma<-P02640_PF02209_Delsub$sigma/sqrt(syns)
  
  #### test if the change in fitness is signifcant
  ## we need to test for significance for values > and < than 1.  
  ## We calculate a z-stat and then do a 2-tailed test
  # we do mulitple testing correction using bonferroni

  # CX
  zstats=(P02640_PF02209_CX$scaled_fitness-1)/P02640_PF02209_CX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P02640_PF02209_CX$zstat=zstats
  P02640_PF02209_CX$pvalue=pvals_man
  P02640_PF02209_CX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P02640_PF02209_CX$fdrsig=P02640_PF02209_CX$FDR<=.05
  P02640_PF02209_CX<-P02640_PF02209_CX %>% mutate(significant=
                                                    ifelse(fdrsig == "TRUE", "*",NA))
  
  
  ##CXX
  zstats=(P02640_PF02209_CXX$scaled_fitness-1)/P02640_PF02209_CXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P02640_PF02209_CXX$zstat=zstats
  P02640_PF02209_CXX$pvalue=pvals_man
  P02640_PF02209_CXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P02640_PF02209_CXX$fdrsig=P02640_PF02209_CXX$FDR<=.05
  P02640_PF02209_CXX<-P02640_PF02209_CXX %>% mutate(significant=
                                                      ifelse(fdrsig == "TRUE", "*",NA))
  ##CXXX
  zstats=(P02640_PF02209_CXXX$scaled_fitness-1)/P02640_PF02209_CXXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P02640_PF02209_CXXX$zstat=zstats
  P02640_PF02209_CXXX$pvalue=pvals_man
  P02640_PF02209_CXXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P02640_PF02209_CXXX$fdrsig=P02640_PF02209_CXXX$FDR<=.05
  P02640_PF02209_CXXX<-P02640_PF02209_CXXX %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #singleDEL
  zstats=(P02640_PF02209_singleDEL$scaled_fitness-1)/P02640_PF02209_singleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P02640_PF02209_singleDEL$zstat=zstats
  P02640_PF02209_singleDEL$pvalue=pvals_man
  P02640_PF02209_singleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P02640_PF02209_singleDEL$fdrsig=P02640_PF02209_singleDEL$FDR<=.05
  P02640_PF02209_singleDEL<-P02640_PF02209_singleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #doubleDEL
  zstats=(P02640_PF02209_doubleDEL$scaled_fitness-1)/P02640_PF02209_doubleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P02640_PF02209_doubleDEL$zstat=zstats
  P02640_PF02209_doubleDEL$pvalue=pvals_man
  P02640_PF02209_doubleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P02640_PF02209_doubleDEL$fdrsig=P02640_PF02209_doubleDEL$FDR<=.05
  P02640_PF02209_doubleDEL<-P02640_PF02209_doubleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #tripleDEL
  zstats=(P02640_PF02209_tripleDEL$scaled_fitness-1)/P02640_PF02209_tripleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P02640_PF02209_tripleDEL$zstat=zstats
  P02640_PF02209_tripleDEL$pvalue=pvals_man
  P02640_PF02209_tripleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P02640_PF02209_tripleDEL$fdrsig=P02640_PF02209_tripleDEL$FDR<=.05
  P02640_PF02209_tripleDEL<-P02640_PF02209_tripleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #DelSub
  zstats=(P02640_PF02209_Delsub$scaled_fitness-1)/P02640_PF02209_Delsub$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  P02640_PF02209_Delsub$zstat=zstats
  P02640_PF02209_Delsub$pvalue=pvals_man
  P02640_PF02209_Delsub[,FDR:=p.adjust(pvalue,method="bonferroni")]
  P02640_PF02209_Delsub$fdrsig=P02640_PF02209_Delsub$FDR<=.05
  P02640_PF02209_Delsub<-P02640_PF02209_Delsub %>% mutate(significant=
                                                            ifelse(fdrsig == "TRUE", "*",NA))
  
  ##some of the DelSubs can cause 2 different substitutions at the same positions. We need to seperate them out into 2 dataframes. 
  P02640_PF02209_Delsub_2<-P02640_PF02209_Delsub[duplicated(P02640_PF02209_Delsub$Pos),]
  P02640_PF02209_Delsub<-P02640_PF02209_Delsub[!duplicated(P02640_PF02209_Delsub$Pos),]
  P02640_PF02209_Delsub_2$type<-"delSub_1"
  P02640_PF02209_Delsub$type<-"delSub_2"
  P02640_PF02209_Delsub_2$Pos<-as.numeric(P02640_PF02209_Delsub_2$Pos)
  P02640_PF02209_Delsub$Pos<-as.numeric(P02640_PF02209_Delsub$Pos)
  
  
  P02640_PF02209_Delsub<-P02640_PF02209_Delsub[P02640_PF02209_Delsub$STOP==FALSE,]
  
  
  
  
  
  ## final df
  P02640_PF02209_allvariants<-rbind(P02640_PF02209_syns[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
                                    P02640_PF02209_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P02640_PF02209_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P02640_PF02209_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P02640_PF02209_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P02640_PF02209_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")], 
                                    P02640_PF02209_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P02640_PF02209_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    P02640_PF02209_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    fill=TRUE)
  
  return(list(P02640_PF02209_allvariants = P02640_PF02209_allvariants, 
              P02640_PF02209_syns = P02640_PF02209_syns[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
              P02640_PF02209_CX = P02640_PF02209_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P02640_PF02209_CXX = P02640_PF02209_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P02640_PF02209_CXXX = P02640_PF02209_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P02640_PF02209_singleDEL = P02640_PF02209_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P02640_PF02209_doubleDEL = P02640_PF02209_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P02640_PF02209_tripleDEL = P02640_PF02209_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P02640_PF02209_Delsub = P02640_PF02209_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              P02640_PF02209_Delsub_2 = P02640_PF02209_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")]))
  
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  