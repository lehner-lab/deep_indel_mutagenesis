normalise_grb2_bind_data <- function() {
  
  ##extract position of mutation, add Pos  and length column
  grb2_bind_synonymous$diff<-as.character(grb2_bind_synonymous$diff)
  grb2_bind_synonymous$diff<-str_first_number(grb2_bind_synonymous$diff)
  grb2_bind_synonymous$diff<-as.numeric(grb2_bind_synonymous$diff)
  grb2_bind_synonymous$Pos<-as.numeric(grb2_bind_synonymous$diff)
  grb2_bind_synonymous$lenght<-str_length(grb2_bind_synonymous$aa_seq)
  
  
  ##extract positions as numeric into the Pos column for all df
  grb2_bind_singleDEL$Pos<-as.numeric(grb2_bind_singleDEL$diff)
  
  grb2_bind_doubleDEL$diff<-as.character(grb2_bind_doubleDEL$diff)
  grb2_bind_doubleDEL$Pos<-str_first_number(grb2_bind_doubleDEL$diff)
  grb2_bind_doubleDEL$Pos<-as.numeric(grb2_bind_doubleDEL$Pos)
  
  grb2_bind_tripleDEL$diff<-as.character(grb2_bind_tripleDEL$diff)
  grb2_bind_tripleDEL$Pos<-str_first_number(grb2_bind_tripleDEL$diff)
  grb2_bind_tripleDEL$Pos<-as.numeric(grb2_bind_tripleDEL$Pos)
  
  grb2_bind_allins$Pos<-as.numeric(grb2_bind_allins$diff)
  grb2_bind_CX$Pos<-as.numeric(grb2_bind_CX$diff)
  
  grb2_bind_CXX$diff<-as.character(grb2_bind_CXX$diff)
  grb2_bind_CXX$Pos<-str_first_number(grb2_bind_CXX$diff)
  grb2_bind_CXX$Pos<-as.numeric(grb2_bind_CXX$Pos)
  
  grb2_bind_CXXX$diff<-as.character(grb2_bind_CXXX$diff)
  grb2_bind_CXXX$Pos<-str_first_number(grb2_bind_CXXX$diff)
  grb2_bind_CXXX$Pos<-as.numeric(grb2_bind_CXXX$Pos)
  
  grb2_bind_Delsub$diff<-as.character(grb2_bind_Delsub$diff)
  grb2_bind_Delsub$Pos<-str_first_number(grb2_bind_Delsub$diff)
  
  ### find the aa position of the synonymous mutation
  grb2_bind_synonymous$Pos <- as.character(floor(grb2_bind_synonymous$diff / 3))
  
  ## add mutation type
  grb2_bind_singleDEL$type<-"singleDEL"
  grb2_bind_doubleDEL$type<-"doubleDEL"
  grb2_bind_tripleDEL$type<-"tripleDEL"
  grb2_bind_synonymous$type<-"synonymous"
  grb2_bind_CX$type<-"singleINS"
  grb2_bind_CXX$type<-"doubleINS"
  grb2_bind_CXXX$type<-"tripleINS"
  grb2_bind_Delsub$type<-"delSub"    
  grb2_bind_allins$type<-"allins"
  
  ## data normalisation
  # normalise the data to mode of dead peak of the bimodal distribution of variants
  
  #bind all data
  test<-rbind(grb2_bind_synonymous,
              grb2_bind_allins,
              grb2_bind_singles,
              grb2_bind_CX,
              grb2_bind_CXX,
              grb2_bind_CNN,
              grb2_bind_insAA,
              grb2_bind_insCC,
              grb2_bind_CXXX,
              grb2_bind_CNNN,
              grb2_bind_singleDEL,
              grb2_bind_doubleDEL, 
              grb2_bind_tripleDEL,
              fill=TRUE)
  
  ##find mode of the lower peak in the bimodal distribution --> referred to as STOPs
  STOPs<-min(mlv(na.omit(test$fitness), method="naive"))
  
  ## substract the mode from all variants. 
  grb2_bind_CX$norm_fitness<-grb2_bind_CX$fitness-STOPs
  grb2_bind_CXX$norm_fitness<-grb2_bind_CXX$fitness-STOPs
  grb2_bind_CXXX$norm_fitness<-grb2_bind_CXXX$fitness-STOPs
  grb2_bind_singleDEL$norm_fitness<-grb2_bind_singleDEL$fitness-STOPs
  grb2_bind_doubleDEL$norm_fitness<-grb2_bind_doubleDEL$fitness-STOPs
  grb2_bind_tripleDEL$norm_fitness<-grb2_bind_tripleDEL$fitness-STOPs
  grb2_bind_Delsub$norm_fitness<-grb2_bind_Delsub$fitness-STOPs
  grb2_bind_synonymous$norm_fitness<-grb2_bind_synonymous$fitness-STOPs
  grb2_bind_insAA$norm_fitness<-grb2_bind_insAA$fitness-STOPs
  grb2_bind_insCC$norm_fitness<-grb2_bind_insAA$fitness-STOPs
  grb2_bind_CNN$norm_fitness<-grb2_bind_CNN$fitness-STOPs
  grb2_bind_CNNN$norm_fitness<-grb2_bind_CNNN$fitness-STOPs
  grb2_bind_allins$norm_fitness<-grb2_bind_allins$fitness-STOPs
  grb2_bind_singles$norm_fitness<-grb2_bind_singles$fitness-STOPs
  
  ## find weighted mean of the synonymous variants
  syns<-sum(grb2_bind_synonymous$norm_fitness/(grb2_bind_synonymous$sigma^2))/sum(1/(grb2_bind_synonymous$sigma^2))
  
  ##normalise by dividing variants with the weighted mean of the synonymous variants
  grb2_bind_synonymous$scaled_fitness<-grb2_bind_synonymous$norm_fitness/syns
  grb2_bind_CX$scaled_fitness<-grb2_bind_CX$norm_fitness/syns
  grb2_bind_CXX$scaled_fitness<-grb2_bind_CXX$norm_fitness/syns
  grb2_bind_CXXX$scaled_fitness<-grb2_bind_CXXX$norm_fitness/syns
  grb2_bind_singleDEL$scaled_fitness<-grb2_bind_singleDEL$norm_fitness/syns
  grb2_bind_doubleDEL$scaled_fitness<-grb2_bind_doubleDEL$norm_fitness/syns
  grb2_bind_tripleDEL$scaled_fitness<-grb2_bind_tripleDEL$norm_fitness/syns
  grb2_bind_Delsub$scaled_fitness<-grb2_bind_Delsub$norm_fitness/syns
  grb2_bind_insAA$scaled_fitness<-grb2_bind_insAA$norm_fitness/syns
  grb2_bind_insCC$scaled_fitness<-grb2_bind_insCC$norm_fitness/syns
  grb2_bind_CNN$scaled_fitness<-grb2_bind_CNN$norm_fitness/syns
  grb2_bind_CNNN$scaled_fitness<-grb2_bind_CNNN$norm_fitness/syns
  grb2_bind_allins$scaled_fitness<-grb2_bind_allins$norm_fitness/syns
  grb2_bind_singles$scaled_fitness<-grb2_bind_singles$norm_fitness/syns
  
  ## scale the error by dividing by the same factor: dividing by the squared root of the synonmous fitness
  grb2_bind_synonymous$scaled_sigma<-grb2_bind_synonymous$sigma/sqrt(syns)
  grb2_bind_CX$scaled_sigma<-grb2_bind_CX$sigma/sqrt(syns)
  grb2_bind_CXX$scaled_sigma<-grb2_bind_CXX$sigma/sqrt(syns)
  grb2_bind_CXXX$scaled_sigma<-grb2_bind_CXXX$sigma/sqrt(syns)
  grb2_bind_singleDEL$scaled_sigma<-grb2_bind_singleDEL$sigma/sqrt(syns)
  grb2_bind_doubleDEL$scaled_sigma<-grb2_bind_doubleDEL$sigma/sqrt(syns)
  grb2_bind_tripleDEL$scaled_sigma<-grb2_bind_tripleDEL$sigma/sqrt(syns)
  grb2_bind_Delsub$scaled_sigma<-grb2_bind_Delsub$sigma/sqrt(syns)
  grb2_bind_insAA$scaled_sigma<-grb2_bind_insAA$sigma/sqrt(syns)
  grb2_bind_insCC$scaled_sigma<-grb2_bind_insCC$sigma/sqrt(syns)
  grb2_bind_CNN$scaled_sigma<-grb2_bind_CNN$sigma/sqrt(syns)
  grb2_bind_CNNN$scaled_sigma<-grb2_bind_CNNN$sigma/sqrt(syns)
  grb2_bind_allins$scaled_sigma<-grb2_bind_allins$sigma/sqrt(syns)
  grb2_bind_singles$scaled_sigma<-grb2_bind_singles$sigma/sqrt(syns)
  
  ##fix mutation positions that are called wrong manually. 
  # this is just an artifact of the function used to call the positions of difference and differing aa sequence
  grb2_bind_doubleDEL[grb2_bind_doubleDEL$diff == "c(20, 22)", ]$Pos<-21
  grb2_bind_doubleDEL[grb2_bind_doubleDEL$diff == "c(35, 37)", ]$Pos<-36
  
  grb2_bind_tripleDEL[grb2_bind_tripleDEL$diff == "c(7, 8, 10)", ]$Pos<-8
  grb2_bind_tripleDEL[grb2_bind_tripleDEL$diff == "c(7, 8, 11)", ]$Pos<-9
  grb2_bind_tripleDEL[grb2_bind_tripleDEL$diff == "c(20, 22, 23)", ]$Pos<-21
  grb2_bind_tripleDEL[grb2_bind_tripleDEL$diff == "c(35, 37, 38)", ]$Pos<-36
  
  ## add missing variants to the df
  # these are variants that have the same aa_seq although the mutation is technically in different positons. 
  #singleDEL:20=21,35=36
  grb2_bind_singleDEL_missing<-grb2_bind_singleDEL[c(40,18),]
  grb2_bind_singleDEL_missing$Pos<-c(21,36)
  grb2_bind_singleDEL<-rbind(grb2_bind_singleDEL,
                             grb2_bind_singleDEL_missing)
  
  #doubleDEL: 7=8, 7=9
  grb2_bind_doubleDEL_missing<-grb2_bind_doubleDEL[c(3,3),]
  grb2_bind_doubleDEL_missing$Pos<-c(8,9)
  grb2_bind_doubleDEL<-rbind(grb2_bind_doubleDEL,
                             grb2_bind_doubleDEL_missing)
  
  #tripleDEL:13=14, 15=16, 29=30, 42=43
  grb2_bind_tripleDEL_missing<-grb2_bind_tripleDEL[c(8, 41, 16, 24),]
  grb2_bind_tripleDEL_missing$Pos<-c(14, 16, 30,43)
  grb2_bind_tripleDEL<-rbind(grb2_bind_tripleDEL,
                             grb2_bind_tripleDEL_missing)
  #CX 20=21, 35=36
  grb2_bind_CX_missing<-grb2_bind_CX[c(11,35),]
  grb2_bind_CX_missing$Pos<-c(21,36)
  grb2_bind_CX<-rbind(grb2_bind_CX,
                      grb2_bind_CX_missing)
  
  #### test if the change in fitness is signifcant
  ## we need to test for significance for values > and < than 1.  
  ## We calculate a z-stat and then do a 2-tailed test
  # we do mulitple testing correction using bonferroni
  
  # CX
  zstats=(grb2_bind_CX$scaled_fitness-1)/grb2_bind_CX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_bind_CX$zstat=zstats
  grb2_bind_CX$pvalue=pvals_man
  grb2_bind_CX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_bind_CX$fdrsig=grb2_bind_CX$FDR<=.05
  grb2_bind_CX<-grb2_bind_CX %>% mutate(significant=
                                          ifelse(fdrsig == "TRUE", "*",NA))
  
  
  ##CXX
  zstats=(grb2_bind_CXX$scaled_fitness-1)/grb2_bind_CXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_bind_CXX$zstat=zstats
  grb2_bind_CXX$pvalue=pvals_man
  grb2_bind_CXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_bind_CXX$fdrsig=grb2_bind_CXX$FDR<=.05
  grb2_bind_CXX<-grb2_bind_CXX %>% mutate(significant=
                                            ifelse(fdrsig == "TRUE", "*",NA))
  ##CXXX
  zstats=(grb2_bind_CXXX$scaled_fitness-1)/grb2_bind_CXXX$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_bind_CXXX$zstat=zstats
  grb2_bind_CXXX$pvalue=pvals_man
  grb2_bind_CXXX[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_bind_CXXX$fdrsig=grb2_bind_CXXX$FDR<=.05
  grb2_bind_CXXX<-grb2_bind_CXXX %>% mutate(significant=
                                              ifelse(fdrsig == "TRUE", "*",NA))
  #singleDEL
  zstats=(grb2_bind_singleDEL$scaled_fitness-1)/grb2_bind_singleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_bind_singleDEL$zstat=zstats
  grb2_bind_singleDEL$pvalue=pvals_man
  grb2_bind_singleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_bind_singleDEL$fdrsig=grb2_bind_singleDEL$FDR<=.05
  grb2_bind_singleDEL<-grb2_bind_singleDEL %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #doubleDEL
  zstats=(grb2_bind_doubleDEL$scaled_fitness-1)/grb2_bind_doubleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_bind_doubleDEL$zstat=zstats
  grb2_bind_doubleDEL$pvalue=pvals_man
  grb2_bind_doubleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_bind_doubleDEL$fdrsig=grb2_bind_doubleDEL$FDR<=.05
  grb2_bind_doubleDEL<-grb2_bind_doubleDEL %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #tripleDEL
  zstats=(grb2_bind_tripleDEL$scaled_fitness-1)/grb2_bind_tripleDEL$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_bind_tripleDEL$zstat=zstats
  grb2_bind_tripleDEL$pvalue=pvals_man
  grb2_bind_tripleDEL[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_bind_tripleDEL$fdrsig=grb2_bind_tripleDEL$FDR<=.05
  grb2_bind_tripleDEL<-grb2_bind_tripleDEL %>% mutate(significant=
                                                        ifelse(fdrsig == "TRUE", "*",NA))
  #DelSub
  zstats=(grb2_bind_Delsub$scaled_fitness-1)/grb2_bind_Delsub$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_bind_Delsub$zstat=zstats
  grb2_bind_Delsub$pvalue=pvals_man
  grb2_bind_Delsub[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_bind_Delsub$fdrsig=grb2_bind_Delsub$FDR<=.05
  grb2_bind_Delsub<-grb2_bind_Delsub %>% mutate(significant=
                                                  ifelse(fdrsig == "TRUE", "*",NA))
  
  grb2_bind_Delsub<-grb2_bind_Delsub[grb2_bind_Delsub$STOP==FALSE,]
  
  
  ##some of the DelSubs can cause 2 different substitutions at the same positions. We need to seperate them out into 2 dataframes. 
  grb2_bind_Delsub_2<-grb2_bind_Delsub[duplicated(grb2_bind_Delsub$Pos),]
  grb2_bind_Delsub<-grb2_bind_Delsub[!duplicated(grb2_bind_Delsub$Pos),]
  grb2_bind_Delsub_2$type<-"delSub_1"
  grb2_bind_Delsub$type<-"delSub_2"
  
  #### CNN
  zstats=(grb2_bind_CNN$scaled_fitness-1)/grb2_bind_CNN$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_bind_CNN$zstat=zstats
  grb2_bind_CNN$pvalue=pvals_man
  grb2_bind_CNN[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_bind_CNN$fdrsig=grb2_bind_CNN$FDR<=.05
  grb2_bind_CNN<-grb2_bind_CNN %>% mutate(significant=
                                            ifelse(fdrsig == "TRUE", "*",NA))
  
  grb2_bind_CNN$Pos<-str_first_number(grb2_bind_CNN$diff)
  
  #### insAA
  zstats=(grb2_bind_insAA$scaled_fitness-1)/grb2_bind_insAA$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_bind_insAA$zstat=zstats
  grb2_bind_insAA$pvalue=pvals_man
  grb2_bind_insAA[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_bind_insAA$fdrsig=grb2_bind_insAA$FDR<=.05
  grb2_bind_insAA<-grb2_bind_insAA %>% mutate(significant=
                                                ifelse(fdrsig == "TRUE", "*",NA))
  
  grb2_bind_insAA$Pos<-str_first_number(grb2_bind_insAA$diff)
  grb2_bind_insAA$type<-"doubleAA"
  
  ### insCC
  zstats=(grb2_bind_insCC$scaled_fitness-1)/grb2_bind_insCC$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_bind_insCC$zstat=zstats
  grb2_bind_insCC$pvalue=pvals_man
  grb2_bind_insCC[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_bind_insCC$fdrsig=grb2_bind_insCC$FDR<=.05
  grb2_bind_insCC<-grb2_bind_insCC %>% mutate(significant=
                                                ifelse(fdrsig == "TRUE", "*",NA))
  
  grb2_bind_insCC$Pos<-str_first_number(grb2_bind_insCC$diff)
  grb2_bind_insCC$type<-"doubleCC"
  
  ### CNNN
  zstats=(grb2_bind_CNNN$scaled_fitness-1)/grb2_bind_CNNN$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_bind_CNNN$zstat=zstats
  grb2_bind_CNNN$pvalue=pvals_man
  grb2_bind_CNNN[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_bind_CNNN$fdrsig=grb2_bind_CNNN$FDR<=.05
  grb2_bind_CNNN<-grb2_bind_CNNN %>% mutate(significant=
                                              ifelse(fdrsig == "TRUE", "*",NA))
  
  grb2_bind_CNNN$Pos<-str_first_number(grb2_bind_CNNN$diff)
  
  ### allins
  zstats=(grb2_bind_allins$scaled_fitness-1)/grb2_bind_allins$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_bind_allins$zstat=zstats
  grb2_bind_allins$pvalue=pvals_man
  grb2_bind_allins[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_bind_allins$fdrsig=grb2_bind_allins$FDR<=.05
  grb2_bind_allins<-grb2_bind_allins %>% mutate(significant=
                                                  ifelse(fdrsig == "TRUE", "*",NA))
  ### singles
  zstats=(grb2_bind_singles$scaled_fitness-1)/grb2_bind_singles$scaled_sigma
  pvals_man=sapply(zstats, function(x){
    if(x <0){
      2*pnorm(q = x, lower.tail = T) 
    }else{
      2*pnorm(q = x, lower.tail = F) 
    }
  })
  grb2_bind_singles$zstat=zstats
  grb2_bind_singles$pvalue=pvals_man
  grb2_bind_singles[,FDR:=p.adjust(pvalue,method="bonferroni")]
  grb2_bind_singles$fdrsig=grb2_bind_singles$FDR<=.05
  grb2_bind_singles<-grb2_bind_singles %>% mutate(significant=
                                                    ifelse(fdrsig == "TRUE", "*",NA))
  grb2_bind_singles$type<-"substitution"
  
  ## remove stop codons 
  grb2_bind_singles <- grb2_bind_singles[grb2_bind_singles$STOP == FALSE,]
  
  ## bind all variants. 
  grb2_bind_allvariants<-rbind(grb2_bind_synonymous[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
                               grb2_bind_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               grb2_bind_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               grb2_bind_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               grb2_bind_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               grb2_bind_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")], 
                               grb2_bind_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               grb2_bind_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               grb2_bind_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
                               grb2_bind_singles[grb2_bind_singles$STOP=="FALSE",c("aa_seq", "Pos", "Mut","type", "scaled_fitness", "scaled_sigma", "significant")],
                               fill=TRUE)
  
  
  return(list(grb2_bind_allvariants = grb2_bind_allvariants, 
              grb2_bind_CX = grb2_bind_CX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_bind_CXX = grb2_bind_CXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_bind_CXXX = grb2_bind_CXXX[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_bind_CNN = grb2_bind_CNN[,c("aa_seq", "Pos","insID","lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_bind_CNNN = grb2_bind_CNNN[,c("aa_seq", "Pos", "insID","lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_bind_insCC = grb2_bind_insCC[,c("aa_seq", "Pos", "insID","lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_bind_insAA = grb2_bind_insAA[,c("aa_seq", "Pos", "insID","lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_bind_allins = grb2_bind_allins[,c("aa_seq", "Pos", "insID","lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_bind_synonymous = grb2_bind_synonymous[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma")],
              grb2_bind_singles = grb2_bind_singles[,c("aa_seq", "Pos", "Mut", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_bind_Delsub = grb2_bind_Delsub[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_bind_Delsub_2 = grb2_bind_Delsub_2[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_bind_singleDEL = grb2_bind_singleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_bind_doubleDEL = grb2_bind_doubleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")],
              grb2_bind_tripleDEL = grb2_bind_tripleDEL[,c("aa_seq", "Pos", "lenght", "type", "scaled_fitness", "scaled_sigma", "significant")]))
}



