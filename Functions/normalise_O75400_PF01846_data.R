normalise_O75400_PF01846_data <- function() {
  
  ##isolate only the programmed lenghts and split into different dfs depending on indel type
  O75400_PF01846_variants$lenght=str_length(O75400_PF01846_variants$aa_seq)
  
  O75400_PF01846_singleDEL <- O75400_PF01846_variants[O75400_PF01846_variants$lenght == "49",]
  O75400_PF01846_doubleDEL <- O75400_PF01846_variants[O75400_PF01846_variants$lenght == "48",]
  O75400_PF01846_tripleDEL <- O75400_PF01846_variants[O75400_PF01846_variants$lenght == "47",]
  O75400_PF01846_CX <- O75400_PF01846_variants[O75400_PF01846_variants$lenght == "51",]
  O75400_PF01846_CXX <- O75400_PF01846_variants[O75400_PF01846_variants$lenght == "52",]
  O75400_PF01846_CXXX <- O75400_PF01846_variants[O75400_PF01846_variants$lenght == "53",]
  
  O75400_PF01846_singleDEL<-O75400_PF01846_singleDEL[O75400_PF01846_singleDEL$characters<3,]
  
  ##Delsubs
  O75400_PF01846_Delsub<-O75400_PF01846_variants[O75400_PF01846_variants$lenght == "49",]
  O75400_PF01846_Delsub<-O75400_PF01846_Delsub[O75400_PF01846_Delsub$characters>2,]
  
  ##extract position of mutation, add Pos  and length column
  O75400_PF01846_syns$diff<-as.character(O75400_PF01846_syns$diff)
  O75400_PF01846_syns$diff<-str_first_number(O75400_PF01846_syns$diff)
  O75400_PF01846_syns$diff<-as.numeric(O75400_PF01846_syns$diff)
  O75400_PF01846_syns$Pos<-as.numeric(O75400_PF01846_syns$diff)
  O75400_PF01846_syns$lenght<-str_length(O75400_PF01846_syns$aa_seq)

  ##extract positions as numeric into the Pos column for all df
  O75400_PF01846_singleDEL$Pos<-as.numeric(O75400_PF01846_singleDEL$diff)
  
  O75400_PF01846_doubleDEL$diff<-as.character(O75400_PF01846_doubleDEL$diff)
  O75400_PF01846_doubleDEL$Pos<-str_first_number(O75400_PF01846_doubleDEL$diff)
  O75400_PF01846_doubleDEL$Pos<-as.numeric(O75400_PF01846_doubleDEL$Pos)
  
  O75400_PF01846_tripleDEL$diff<-as.character(O75400_PF01846_tripleDEL$diff)
  O75400_PF01846_tripleDEL$Pos<-str_first_number(O75400_PF01846_tripleDEL$diff)
  O75400_PF01846_tripleDEL$Pos<-as.numeric(O75400_PF01846_tripleDEL$Pos)
  
  O75400_PF01846_CX$Pos<-as.numeric(O75400_PF01846_CX$diff)
  
  O75400_PF01846_CXX$diff<-as.character(O75400_PF01846_CXX$diff)
  O75400_PF01846_CXX$Pos<-str_first_number(O75400_PF01846_CXX$diff)
  O75400_PF01846_CXX$Pos<-as.numeric(O75400_PF01846_CXX$Pos)
  
  O75400_PF01846_CXXX$diff<-as.character(O75400_PF01846_CXXX$diff)
  O75400_PF01846_CXXX$Pos<-str_first_number(O75400_PF01846_CXXX$diff)
  O75400_PF01846_CXXX$Pos<-as.numeric(O75400_PF01846_CXXX$Pos)
  
  O75400_PF01846_Delsub$diff<-as.character(O75400_PF01846_Delsub$diff)
  O75400_PF01846_Delsub$Pos<-str_first_number(O75400_PF01846_Delsub$diff)
  
  ### find the aa position of the synonymous mutation
  O75400_PF01846_syns$Pos <- as.character(floor(O75400_PF01846_syns$diff / 3))
  
  
  ## add mutation type
  O75400_PF01846_singleDEL$type<-"singleDEL"
  O75400_PF01846_doubleDEL$type<-"doubleDEL"
  O75400_PF01846_tripleDEL$type<-"tripleDEL"
  O75400_PF01846_CX$type<-"singleINS"
  O75400_PF01846_CXX$type<-"doubleINS"
  O75400_PF01846_CXXX$type<-"tripleINS"
  O75400_PF01846_Delsub$type<-"delSub"    
  O75400_PF01846_syns$type<-"synonymous"
  
  ##fix mutation positions that are called wrong manually. 
  # this is just an artifact of the function used to call the positions of difference and differing aa sequence
  O75400_PF01846_doubleDEL[O75400_PF01846_doubleDEL$diff=="c(9, 11)",]$Pos<-10
  O75400_PF01846_doubleDEL[O75400_PF01846_doubleDEL$diff=="c(28, 30)",]$Pos<-29
  O75400_PF01846_doubleDEL[O75400_PF01846_doubleDEL$diff=="c(43, 45)",]$Pos<-44
  
  O75400_PF01846_tripleDEL[O75400_PF01846_tripleDEL$diff=="c(9, 11, 12)",]$Pos<-10
  O75400_PF01846_tripleDEL[O75400_PF01846_tripleDEL$diff=="c(11, 12, 14)",]$Pos<-12
  O75400_PF01846_tripleDEL[O75400_PF01846_tripleDEL$diff=="c(25, 26, 28)",]$Pos<-26
  O75400_PF01846_tripleDEL[O75400_PF01846_tripleDEL$diff=="c(28, 30, 31)",]$Pos<-29
  O75400_PF01846_tripleDEL[O75400_PF01846_tripleDEL$diff=="c(43, 45, 46)",]$Pos<-44
  O75400_PF01846_tripleDEL[O75400_PF01846_tripleDEL$diff=="c(36, 37, 39)",]$Pos<-38
  
  O75400_PF01846_CXX[O75400_PF01846_CXX$diff=="c(9, 11)",]$Pos<-10
  O75400_PF01846_CXX[O75400_PF01846_CXX$diff=="c(28, 30)",]$Pos<-29
  O75400_PF01846_CXX[O75400_PF01846_CXX$diff=="c(43, 45)",]$Pos<-44
  
  O75400_PF01846_CXXX[O75400_PF01846_CXXX$diff=="c(9, 11, 12)",]$Pos<-10
  O75400_PF01846_CXXX[O75400_PF01846_CXXX$diff=="c(11, 12, 14)",]$Pos<-12
  O75400_PF01846_CXXX[O75400_PF01846_CXXX$diff=="c(25, 26, 28)",]$Pos<-26
  O75400_PF01846_CXXX[O75400_PF01846_CXXX$diff=="c(28, 30, 31)",]$Pos<-29
  O75400_PF01846_CXXX[O75400_PF01846_CXXX$diff=="c(43, 45, 46)",]$Pos<-44
  O75400_PF01846_CXXX[O75400_PF01846_CXXX$diff=="c(36, 37, 39)",]$Pos<-38
  
  ## add missing variants to the df
  # these are variants that have the same aa_seq although the mutation is technically in different positons. 
  #singleDEL: 9=10; 28=29, 43=44, 
  O75400_PF01846_singleDEL_missing<-O75400_PF01846_singleDEL[c(3,11,26),]
  O75400_PF01846_singleDEL_missing$Pos<-c(10, 29, 44)
  O75400_PF01846_singleDEL<-rbind(O75400_PF01846_singleDEL,
                                  O75400_PF01846_singleDEL_missing)
  #doubleDEL: 11=12; 25=26, 36=37)
  O75400_PF01846_doubleDEL_missing<-O75400_PF01846_doubleDEL[c(7, 36, 17),]
  O75400_PF01846_doubleDEL_missing$Pos<-c(12,26,37)
  O75400_PF01846_doubleDEL<-rbind(O75400_PF01846_doubleDEL,
                                  O75400_PF01846_doubleDEL_missing)
  #tripleDEL: 2=3,17=18, 36=37, 46=47
  O75400_PF01846_tripleDEL_missing<-O75400_PF01846_tripleDEL[c(44,37,17,22),]
  O75400_PF01846_tripleDEL_missing$Pos<-c(3,18,37,47)
  O75400_PF01846_tripleDEL<-rbind(O75400_PF01846_tripleDEL,
                                  O75400_PF01846_tripleDEL_missing)
  #singleCX:  9=10; 28=29, 43=44, 
  O75400_PF01846_CX_missing<-O75400_PF01846_CX[c(45, 37, 22),]
  O75400_PF01846_CX_missing$Pos<-c(10, 29, 44)
  O75400_PF01846_CX<-rbind(O75400_PF01846_CX,
                           O75400_PF01846_CX_missing)
  #doubleCXX:  11=12; 25=26, 36=37)
  O75400_PF01846_CXX_missing<-O75400_PF01846_CXX[c(40,11,30),]
  O75400_PF01846_CXX_missing$Pos<-c(12, 26,37)
  O75400_PF01846_CXX<-rbind(O75400_PF01846_CXX,
                            O75400_PF01846_CXX_missing)
  #tripleCXX: 2=3; 17=18, 37=38, 46=47  
  O75400_PF01846_CXXX_missing<-O75400_PF01846_CXXX[c(1,8,18,23),]
  O75400_PF01846_CXXX_missing$Pos<-c(3,18,37,47)
  O75400_PF01846_CXXX<-rbind(O75400_PF01846_CXXX,
                             O75400_PF01846_CXXX_missing)
  
  ## data normalisation
  # normalise the data to mode of dead peak of the bimodal distribution of variants
  
  #bind all data
  test<-rbind(O75400_PF01846_syns,
              O75400_PF01846_CX,
              O75400_PF01846_CXX,
              O75400_PF01846_CXXX,
              O75400_PF01846_singleDEL,
              O75400_PF01846_doubleDEL, 
              O75400_PF01846_tripleDEL,
              fill=TRUE)
  
  ##find mode of the lower peak in the bimodal distribution --> referred to as STOPs
  STOPs<-mlv(na.omit(test$fitness), method="naive")
  
  ## substract the mode from all variants. 
  O75400_PF01846_CX$norm_fitness<-O75400_PF01846_CX$fitness-STOPs
  O75400_PF01846_CXX$norm_fitness<-O75400_PF01846_CXX$fitness-STOPs
  O75400_PF01846_CXXX$norm_fitness<-O75400_PF01846_CXXX$fitness-STOPs
  O75400_PF01846_singleDEL$norm_fitness<-O75400_PF01846_singleDEL$fitness-STOPs
  O75400_PF01846_doubleDEL$norm_fitness<-O75400_PF01846_doubleDEL$fitness-STOPs
  O75400_PF01846_tripleDEL$norm_fitness<-O75400_PF01846_tripleDEL$fitness-STOPs
  O75400_PF01846_Delsub$norm_fitness<-O75400_PF01846_Delsub$fitness-STOPs
  O75400_PF01846_syns$norm_fitness<-O75400_PF01846_syns$fitness-STOPs
  
  ## find weighted mean of the synonymous variants
  syns<-sum(O75400_PF01846_syns$norm_fitness/(O75400_PF01846_syns$sigma^2))/sum(1/(O75400_PF01846_syns$sigma^2))
  
  ##normalise by dividing variants with the weighted mean of the synonymous variants
  O75400_PF01846_syns$scaled_fitness<-O75400_PF01846_syns$norm_fitness/syns
  O75400_PF01846_CX$scaled_fitness<-O75400_PF01846_CX$norm_fitness/syns
  O75400_PF01846_CXX$scaled_fitness<-O75400_PF01846_CXX$norm_fitness/syns
  O75400_PF01846_CXXX$scaled_fitness<-O75400_PF01846_CXXX$norm_fitness/syns
  O75400_PF01846_singleDEL$scaled_fitness<-O75400_PF01846_singleDEL$norm_fitness/syns
  O75400_PF01846_doubleDEL$scaled_fitness<-O75400_PF01846_doubleDEL$norm_fitness/syns
  O75400_PF01846_tripleDEL$scaled_fitness<-O75400_PF01846_tripleDEL$norm_fitness/syns
  O75400_PF01846_Delsub$scaled_fitness<-O75400_PF01846_Delsub$norm_fitness/syns
  
  
  ## scale the error by dividing by the same factor: dividing by the squared root of the synonmous fitness
  O75400_PF01846_syns$scaled_sigma<-O75400_PF01846_syns$sigma/sqrt(syns)
  O75400_PF01846_CX$scaled_sigma<-O75400_PF01846_CX$sigma/sqrt(syns)
  O75400_PF01846_CXX$scaled_sigma<-O75400_PF01846_CXX$sigma/sqrt(syns)
  O75400_PF01846_CXXX$scaled_sigma<-O75400_PF01846_CXXX$sigma/sqrt(syns)
  O75400_PF01846_singleDEL$scaled_sigma<-O75400_PF01846_singleDEL$sigma/sqrt(syns)
  O75400_PF01846_doubleDEL$scaled_sigma<-O75400_PF01846_doubleDEL$sigma/sqrt(syns)
  O75400_PF01846_tripleDEL$scaled_sigma<-O75400_PF01846_tripleDEL$sigma/sqrt(syns)
  O75400_PF01846_Delsub$scaled_sigma<-O75400_PF01846_Delsub$sigma/sqrt(syns)
  
  #### test if the change in fitness is signifcant
  ## we need to test for significance for values > and < than 1.  
  ## We calculate a z-stat and then do a 2-tailed test
  # we do mulitple testing correction using bonferroni
  
  # CX
  zstats=(O75400_PF01846_CX$scaled_fitness-1)/O75400_PF01846_CX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  
  O75400_PF01846_CX$zstat=zstats
  O75400_PF01846_CX$pvalue=pvals_man
  O75400_PF01846_CX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  O75400_PF01846_CX$fdrsig=O75400_PF01846_CX$FDR<=.05
  O75400_PF01846_CX<-O75400_PF01846_CX %>% mutate(significant=
                                                    ifelse(fdrsig == "TRUE", "*",NA))
  
  
  ##CXX
  zstats=(O75400_PF01846_CXX$scaled_fitness-1)/O75400_PF01846_CXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  O75400_PF01846_CXX$zstat=zstats
  O75400_PF01846_CXX$pvalue=pvals_man
  O75400_PF01846_CXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  O75400_PF01846_CXX$fdrsig=O75400_PF01846_CXX$FDR<=.05
  O75400_PF01846_CXX<-O75400_PF01846_CXX %>% mutate(significant=
                                                      ifelse(fdrsig == "TRUE", "*",NA))
  ##CXXX
  zstats=(O75400_PF01846_CXXX$scaled_fitness-1)/O75400_PF01846_CXXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  O75400_PF01846_CXXX$zstat=zstats
  O75400_PF01846_CXXX$pvalue=pvals_man
  O75400_PF01846_CXXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  O75400_PF01846_CXXX$fdrsig=O75400_PF01846_CXXX$FDR<=.05
  O75400_PF01846_CXXX<-O75400_PF01846_CXXX %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #singleDEL
  zstats=(O75400_PF01846_singleDEL$scaled_fitness-1)/O75400_PF01846_singleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  O75400_PF01846_singleDEL$zstat=zstats
  O75400_PF01846_singleDEL$pvalue=pvals_man
  O75400_PF01846_singleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  O75400_PF01846_singleDEL$fdrsig=O75400_PF01846_singleDEL$FDR<=.05
  O75400_PF01846_singleDEL<-O75400_PF01846_singleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #doubleDEL
  zstats=(O75400_PF01846_doubleDEL$scaled_fitness-1)/O75400_PF01846_doubleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  O75400_PF01846_doubleDEL$zstat=zstats
  O75400_PF01846_doubleDEL$pvalue=pvals_man
  O75400_PF01846_doubleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  O75400_PF01846_doubleDEL$fdrsig=O75400_PF01846_doubleDEL$FDR<=.05
  O75400_PF01846_doubleDEL<-O75400_PF01846_doubleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #tripleDEL
  zstats=(O75400_PF01846_tripleDEL$scaled_fitness-1)/O75400_PF01846_tripleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  O75400_PF01846_tripleDEL$zstat=zstats
  O75400_PF01846_tripleDEL$pvalue=pvals_man
  O75400_PF01846_tripleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  O75400_PF01846_tripleDEL$fdrsig=O75400_PF01846_tripleDEL$FDR<=.05
  O75400_PF01846_tripleDEL<-O75400_PF01846_tripleDEL %>% mutate(significant=
                                                                  ifelse(fdrsig == "TRUE", "*",NA))
  #DelSub
  zstats=(O75400_PF01846_Delsub$scaled_fitness-1)/O75400_PF01846_Delsub$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  O75400_PF01846_Delsub$zstat=zstats
  O75400_PF01846_Delsub$pvalue=pvals_man
  O75400_PF01846_Delsub[,FDR:=p.adjust(pvalue,method="bonferroni")]
  O75400_PF01846_Delsub$fdrsig=O75400_PF01846_Delsub$FDR<=.05
  O75400_PF01846_Delsub<-O75400_PF01846_Delsub %>% mutate(significant=
                                                            ifelse(fdrsig == "TRUE", "*",NA))
  
  ##some of the DelSubs can cause 2 different substitutions at the same positions. We need to seperate them out into 2 dataframes. 
  O75400_PF01846_Delsub_2<-O75400_PF01846_Delsub[duplicated(O75400_PF01846_Delsub$Pos),]
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  O75400_PF01846_Delsub<-O75400_PF01846_Delsub[!duplicated(O75400_PF01846_Delsub$Pos),]
  O75400_PF01846_Delsub_2$type<-"delSub_1"
  O75400_PF01846_Delsub$type<-"delSub_2"
  O75400_PF01846_Delsub_2$Pos<-as.numeric(O75400_PF01846_Delsub_2$Pos)
  O75400_PF01846_Delsub$Pos<-as.numeric(O75400_PF01846_Delsub$Pos)
  
  
  O75400_PF01846_Delsub<-O75400_PF01846_Delsub[O75400_PF01846_Delsub$STOP==FALSE,]
  
  
  
  
  ## final df
  O75400_PF01846_allvariants<-rbind(O75400_PF01846_syns[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
                                    O75400_PF01846_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    O75400_PF01846_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    O75400_PF01846_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    O75400_PF01846_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    O75400_PF01846_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")], 
                                    O75400_PF01846_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    O75400_PF01846_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    O75400_PF01846_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                                    fill=TRUE)
  
  return(list(O75400_PF01846_allvariants = O75400_PF01846_allvariants, 
              O75400_PF01846_syns = O75400_PF01846_syns[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
              O75400_PF01846_CX = O75400_PF01846_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              O75400_PF01846_CXX = O75400_PF01846_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              O75400_PF01846_CXXX = O75400_PF01846_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              O75400_PF01846_singleDEL = O75400_PF01846_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              O75400_PF01846_doubleDEL = O75400_PF01846_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              O75400_PF01846_tripleDEL = O75400_PF01846_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              O75400_PF01846_Delsub = O75400_PF01846_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              O75400_PF01846_Delsub_2 = O75400_PF01846_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")]))
  
}


  
  
  
  
  