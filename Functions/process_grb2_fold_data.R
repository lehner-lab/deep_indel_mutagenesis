### this script defines the function to call the different mutational variants for grb2-sh3
process_grb2_fold_data <- function() {
  
  load(paste(dimsum_file_path, "grb2_fold_fitness_replicates.RData",sep=""))

  grb2_fold_all_variants<-all_variants
  grb2_fold_singles<-singles
  grb2_fold_synonymous<-synonymous
  grb2_fold_wildtype<-wildtype
  
  # Identify positions where sequences differ from wt_seq at the aa level
  grb2_fold_all_variants$diff=mapply(function(x, y) 
    grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "TYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYV",
    grb2_fold_all_variants$aa_seq)
  
  # Add length of the sequence to the data frame
  grb2_fold_all_variants$lenght=str_length(grb2_fold_all_variants$aa_seq)
  
  # Count the number of characters where sequences differ from wt_seq
  grb2_fold_all_variants$characters<-nchar(grb2_fold_all_variants$diff)
  
  # Isolate the correct sequences
  grb2_fold_variants<-grb2_fold_all_variants[grb2_fold_all_variants$characters <14,]
  grb2_fold_variants<-grb2_fold_variants[grb2_fold_variants$indel==TRUE,]
  
  # Identify positions where sequences differ from wt_seq at the nt level
  grb2_fold_synonymous$diff=mapply(function(x, y) 
    grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "acatacgtccaggccctctttgactttgatccccaggaggatggagagctgggcttccgccggggagattttatccatgtcatggataactcagaccccaactggtggaaaggagcttgccacgggcagaccggcatgtttccccgcaattatgtc",
    grb2_fold_synonymous$nt_seq)
  
  # Count the number of characters where sequences differ from wt_seq
  grb2_fold_synonymous$characters<-nchar(grb2_fold_synonymous$diff)
  
  # Filter the data frame for synonymous mutants
  grb2_fold_synonymous<-grb2_fold_synonymous[grb2_fold_synonymous$characters <12,]
  
  ##isolate only the programmed lenghts and split into different dfs depending on indel type
  grb2_fold_singleDEL <- grb2_fold_variants[grb2_fold_variants$lenght == "51",]
  grb2_fold_doubleDEL <- grb2_fold_variants[grb2_fold_variants$lenght == "50",]
  grb2_fold_tripleDEL <- grb2_fold_variants[grb2_fold_variants$lenght == "49",]
  grb2_fold_CX <- grb2_fold_variants[grb2_fold_variants$lenght == "53",]
  grb2_fold_CXX <- grb2_fold_variants[grb2_fold_variants$lenght == "54",]
  grb2_fold_CXXX <- grb2_fold_variants[grb2_fold_variants$lenght == "55",]
  
  grb2_fold_singleDEL<-grb2_fold_singleDEL[grb2_fold_singleDEL$characters<3,]
  
  grb2_fold_Delsub<-grb2_fold_variants[grb2_fold_variants$lenght == "51",]
  grb2_fold_Delsub<-grb2_fold_Delsub[grb2_fold_Delsub$characters>2,]
  
  ## For the singleINS I need to  extract the single CCC insertions only
  grb2_fold_CX$diff<-as.numeric(grb2_fold_CX$diff)
  grb2_fold_CX$wt_seq<-"TYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYV"
  
  ## Adding the wt_aa column 
  grb2_fold_CX$WT_AA<-substr(grb2_fold_CX$wt_seq, grb2_fold_CX$diff, grb2_fold_CX$diff)
  grb2_fold_CX<-grb2_fold_CX[,-36]
  
  ## Identifying what is the single CCC insertion (refered to as CX) versus all possible insertions
  grb2_fold_CX$insID<-substr(grb2_fold_CX$aa_seq, grb2_fold_CX$diff, grb2_fold_CX$diff)
  grb2_fold_CX$CX<-grb2_fold_CX$WT_AA == grb2_fold_CX$insID
  
  ## Seperating all 20 insertions into a seperate df: grb2_fold_allins and single CCC insertions into: grb2_fold_CX
  grb2_fold_allins<-grb2_fold_CX
  grb2_fold_CX<-grb2_fold_CX[grb2_fold_CX$CX == TRUE,]
  
  ###for the double insertions, I have the double CCC (CXX), the double repeats of the wt_aa (CNN), the double Ala insertions and double Cys insertion 
  ## identify the inserted aa
  grb2_fold_CXX$diff<-as.character(grb2_fold_CXX$diff)
  grb2_fold_CXX$insID<-substr(grb2_fold_CXX$aa_seq, str_first_number(grb2_fold_CXX$diff), str_last_number(grb2_fold_CXX$diff))
  
  ## Split double AA and CC from the rest
  grb2_fold_insAA<-grb2_fold_CXX[grb2_fold_CXX$insID =="AA",]
  grb2_fold_insCC<-grb2_fold_CXX[grb2_fold_CXX$insID =="CC",]
  
  ## Remove double Ala and double Cys from CXX
  grb2_fold_CXX<-grb2_fold_CXX[!grepl("AA", grb2_fold_CXX$insID),]
  grb2_fold_CXX<-grb2_fold_CXX[!grepl("CC", grb2_fold_CXX$insID),]
  
  ## Now split between CNN and CXX by looking for repeated characters in string insID
  grb2_fold_CXX$duplicated<-substr(grb2_fold_CXX$insID, 1, 1) == substr(grb2_fold_CXX$insID, 2, 2)
  
  ## Manually change the variants that are to be present in both CXX and CNN
  missing_from_CXX<-grb2_fold_CXX[c(20,64),]
  
  grb2_fold_CXX[21,]$insID<-"RG"
  grb2_fold_CXX[21,]$duplicated <- "FALSE"
  grb2_fold_CXX[63,]$insID<-"WK"
  grb2_fold_CXX[63,]$duplicated <- "FALSE"
  
  ## seperate CNN and CXX dataframes
  grb2_fold_CNN<-grb2_fold_CXX[grb2_fold_CXX$duplicated == TRUE,]
  grb2_fold_CXX<-grb2_fold_CXX[!grb2_fold_CXX$duplicated == TRUE,]
  
  ## add the missing variants
  grb2_fold_CXX<-rbind(grb2_fold_CXX,
                       missing_from_CXX)
  
  ##add the 9:10 insertion which is same as 7:8 insertion (same aa_seq)
  missing<-grb2_fold_CXX[grb2_fold_CXX$diff == "7:8",]
  missing$diff<-"9:10"
  grb2_fold_CXX<-rbind(grb2_fold_CXX, 
                       missing)
  
  ###re-organise the df based on position
  grb2_fold_CXX <- arrange(grb2_fold_CXX, as.numeric(sub(":.*", "", diff)))
  
  ## manually change to the correct positions
  grb2_fold_CXX[49,]$diff<-"21:22"
  grb2_fold_CXX[50,]$diff<-"36:37"
  
  ###re-organise the df based on position
  grb2_fold_CXX <- arrange(grb2_fold_CXX, as.numeric(sub(":.*", "", diff)))
  
  ### add the 8:9 insertion which is same as 7:8 insertion (same aa_seq)
  missing<-grb2_fold_CXX[c(7),]
  
  missing$diff<-"8:9"
  grb2_fold_CXX<-rbind(grb2_fold_CXX,
                       missing)
  
  ###re-organise the df based on position
  grb2_fold_CXX <- arrange(grb2_fold_CXX, as.numeric(sub(":.*", "", diff)))
  
  ## manually change into the correct insertions
  grb2_fold_CXX[c(7:9),]$insID<-c("FD", "DF", "FD")
  
  ###re-organise the df based on position
  grb2_fold_CNN <- arrange(grb2_fold_CNN, as.numeric(sub(":.*", "", diff)))
  
  ## add missing insertions 20:21 is the same as 21:22 and 35:36 is the same as 36:37 sequence-wise
  missing<-grb2_fold_CNN[c(19,33),]
  missing$diff<-c("21:22", "36:37")
  
  grb2_fold_CNN<-rbind(grb2_fold_CNN,
                       missing)
  
  ## add the missing CNN repeats that are overlapping with the 2 Ala repeats
  missing<-rbind(grb2_fold_insAA[48,],
                 grb2_fold_insAA[14,])
  
  grb2_fold_CNN<-rbind(grb2_fold_CNN,
                       missing, fill=TRUE)
  
  ###re-organise the df based on position
  grb2_fold_CNN <- arrange(grb2_fold_CNN, as.numeric(sub(":.*", "", diff)))
  
  grb2_fold_CNN$type<-"doubleINS_CNN"
  grb2_fold_CXX$type<-"doubleINS_CXX"
  
  
  ### Same is now done for the triple insertions 
  grb2_fold_CXXX$diff<-as.character(grb2_fold_CXXX$diff)
  grb2_fold_CXXX$insID<-substr(grb2_fold_CXXX$aa_seq, str_first_number(grb2_fold_CXXX$diff), str_last_number(grb2_fold_CXXX$diff))
  
  ########################
  ##now split between CNNN and CXXX by looking for repeated characters in string insID
  grb2_fold_CXXX$duplicated<- substr(grb2_fold_CXXX$insID, 1, 2) == substr(grb2_fold_CXXX$insID, 2, 3)
  
  grb2_fold_CNNN<-grb2_fold_CXXX[grb2_fold_CXXX$duplicated == TRUE,]
  grb2_fold_CXXX<-grb2_fold_CXXX[!grb2_fold_CXXX$duplicated == TRUE,]
  
  ##add the missing ones to CXXX
  grb2_fold_CXXX <- arrange(grb2_fold_CXXX, as.numeric(sub(":.*", "", diff)))
  
  missing<-grb2_fold_CXXX[c(43:46),]
  missing$diff<-c("21:23", "36:38", "9:11", "8:10")
  missing$insID<-c("RGD", "WKG", "FDP", "DFD")
  
  grb2_fold_CXXX<-grb2_fold_CXXX[-c(43:46),]
  
  grb2_fold_CXXX<-rbind(grb2_fold_CXXX, 
                        missing)
  
  grb2_fold_CXXX <- arrange(grb2_fold_CXXX, as.numeric(sub(":.*", "", diff)))
  
  missing<-grb2_fold_CXXX[c(13:14),]
  missing$diff<-c("14:15","16:17")
  missing$insID<-c("DGE", "ELG")
  
  grb2_fold_CXXX<-rbind(grb2_fold_CXXX,
                        missing)
  
  grb2_fold_CXXX <- arrange(grb2_fold_CXXX, as.numeric(sub(":.*", "", diff)))
  
  
  missing<-grb2_fold_CXXX[c(29),]
  missing$diff<-c("30:31")
  missing$insID<-c("NSD")
  
  grb2_fold_CXXX<-rbind(grb2_fold_CXXX,
                        missing)
  
  grb2_fold_CXXX <- arrange(grb2_fold_CXXX, as.numeric(sub(":.*", "", diff)))
  
  missing<-grb2_fold_CXXX[c(42),]
  missing$diff<-c("43:44")
  missing$insID<-c("QTG")
  
  grb2_fold_CXXX<-rbind(grb2_fold_CXXX,
                        missing)
  
  grb2_fold_CXXX <- arrange(grb2_fold_CXXX, as.numeric(sub(":.*", "", diff)))
  
  
  ##add the missing ones to CNNN df
  grb2_fold_CNNN <- arrange(grb2_fold_CNNN, as.numeric(sub(":.*", "", diff)))
  
  missing<-grb2_fold_CNNN[c(19,32),]
  missing$diff<-c("21:23", "36:38")
  
  grb2_fold_CNNN<-rbind(grb2_fold_CNNN,
                        missing)
  
  grb2_fold_CNNN <- arrange(grb2_fold_CNNN, as.numeric(sub(":.*", "", diff)))
  
  
  grb2_fold_CNNN$type<-"tripleINS_CNNN"
  grb2_fold_CXXX$type<-"tripleINS_CXXX"
  
  # Return the results
  return(list(grb2_fold_CX = grb2_fold_CX,
              grb2_fold_CXX = grb2_fold_CXX,
              grb2_fold_CNN = grb2_fold_CNN,
              grb2_fold_insAA = grb2_fold_insAA,
              grb2_fold_insCC = grb2_fold_insCC,
              grb2_fold_CXXX = grb2_fold_CXXX,
              grb2_fold_CNNN = grb2_fold_CNNN,
              grb2_fold_singles = grb2_fold_singles,
              grb2_fold_allins = grb2_fold_allins,
              grb2_fold_Delsub = grb2_fold_Delsub,
              grb2_fold_wildtype = grb2_fold_wildtype,
              grb2_fold_synonymous = grb2_fold_synonymous,
              grb2_fold_singleDEL = grb2_fold_singleDEL,
              grb2_fold_doubleDEL = grb2_fold_doubleDEL,
              grb2_fold_tripleDEL = grb2_fold_tripleDEL
              ))
  
}
