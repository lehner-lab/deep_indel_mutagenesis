process_pdz3_fold_data <- function() {
  
  
  load(paste(dimsum_file_path, "pdz3_fold_fitness_replicates.RData",sep=""))

  pdz3_fold_all_variants<-all_variants
  pdz3_fold_singles<-singles
  pdz3_fold_synonymous<-synonymous
  pdz3_fold_wildtype<-wildtype
  
  # Identify positions where sequences differ from wt_seq at the aa level
  pdz3_fold_all_variants$diff=mapply(function(x, y) 
    grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "PRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSV",
    pdz3_fold_all_variants$aa_seq)
  
  # Add length of the sequence to the data frame
  pdz3_fold_all_variants$lenght=str_length(pdz3_fold_all_variants$aa_seq)
  
  # Count the number of characters where sequences differ from wt_seq
  pdz3_fold_all_variants$characters<-nchar(pdz3_fold_all_variants$diff)
  
  # Isolate the correct sequences
  pdz3_variants<-pdz3_fold_all_variants[pdz3_fold_all_variants$characters <14,]
  pdz3_variants<-pdz3_variants[pdz3_variants$indel==TRUE,]
  
  # Identify positions where sequences differ from wt_seq at the nt level
  pdz3_fold_synonymous$diff=mapply(function(x, y) 
    grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "ccgaggcgaattgtgatccaccggggctccacgggcctgggcttcaacatcgtgggtggcgaggacggtgaaggcatcttcatctcctttatcctggccgggggccctgcagacctcagtggggagctgcggaagggggaccagatcctgtcggtc",
    pdz3_fold_synonymous$nt_seq)
  
  # Count the number of characters where sequences differ from wt_seq
  pdz3_fold_synonymous$characters<-nchar(pdz3_fold_synonymous$diff)
  
  # Filter the data frame for synonymous mutants
  pdz3_fold_synonymous<-pdz3_fold_synonymous[pdz3_fold_synonymous$characters <12,]
  
  ##isolate only the programmed lenghts and split into different dfs depending on indel type
  pdz3_fold_singleDEL <- pdz3_variants[pdz3_variants$lenght == "51",]
  pdz3_fold_doubleDEL <- pdz3_variants[pdz3_variants$lenght == "50",]
  pdz3_fold_tripleDEL <- pdz3_variants[pdz3_variants$lenght == "49",]
  pdz3_fold_CX <- pdz3_variants[pdz3_variants$lenght == "53",]
  pdz3_fold_CXX <- pdz3_variants[pdz3_variants$lenght == "54",]
  pdz3_fold_CXXX <- pdz3_variants[pdz3_variants$lenght == "55",]
  
  pdz3_fold_singleDEL<-pdz3_fold_singleDEL[pdz3_fold_singleDEL$characters<3,]
  
  pdz3_fold_Delsub<-pdz3_variants[pdz3_variants$lenght == "51",]
  pdz3_fold_Delsub<-pdz3_fold_Delsub[pdz3_fold_Delsub$characters>2,]
  
  ## For the singleINS I need to  extract the single CCC insertions only
  pdz3_fold_CX$diff<-as.numeric(pdz3_fold_CX$diff)
  pdz3_fold_CX$wt_seq<-"PRRIVIHRGSTGLGFNIVGGEDGEGIFISFILAGGPADLSGELRKGDQILSV"
  
  ## Adding the wt_aa column 
  pdz3_fold_CX$WT_AA<-substr(pdz3_fold_CX$wt_seq, pdz3_fold_CX$diff, pdz3_fold_CX$diff)
  pdz3_fold_CX<-pdz3_fold_CX[,-36]
  
  ## Identifying what is the single CCC insertion (refered to as CX) versus all possible insertions
  pdz3_fold_CX$insID<-substr(pdz3_fold_CX$aa_seq, pdz3_fold_CX$diff, pdz3_fold_CX$diff)
  pdz3_fold_CX$CX<-pdz3_fold_CX$WT_AA == pdz3_fold_CX$insID
  
  ## Seperating all 20 insertions into a seperate df: pdz3_fold_allins and single CCC insertions into: pdz3_fold_CX
  pdz3_fold_allins<-pdz3_fold_CX
  pdz3_fold_CX<-pdz3_fold_CX[pdz3_fold_CX$CX == TRUE,]
  
  
  ###for the double insertions, I have the double CCC (CXX), the double repeats of the wt_aa (CNN), the double Ala insertions and double Cys insertion 
  ## identify the inserted aa
  pdz3_fold_CXX$diff<-as.character(pdz3_fold_CXX$diff)
  pdz3_fold_CXX$insID<-substr(pdz3_fold_CXX$aa_seq, str_first_number(pdz3_fold_CXX$diff), str_last_number(pdz3_fold_CXX$diff))
  
  ## Split double AA and CC from the rest
  pdz3_fold_insAA<-pdz3_fold_CXX[pdz3_fold_CXX$insID =="AA",]
  pdz3_fold_insCC<-pdz3_fold_CXX[pdz3_fold_CXX$insID =="CC",]
  
  ## Remove double Ala and double Cys from CXX
  pdz3_fold_CXX<-pdz3_fold_CXX[!grepl("AA", pdz3_fold_CXX$insID),]
  pdz3_fold_CXX<-pdz3_fold_CXX[!grepl("CC", pdz3_fold_CXX$insID),]
  
  ## Now split between CNN and CXX by looking for repeated characters in string insID
  pdz3_fold_CXX$duplicated<-substr(pdz3_fold_CXX$insID, 1, 1) == substr(pdz3_fold_CXX$insID, 2, 2)

  ## Manually change the variants that are to be present in both CXX and CNN
  pdz3_fold_CNN<-pdz3_fold_CXX[pdz3_fold_CXX$duplicated == TRUE,]
  pdz3_fold_CXX <- arrange(pdz3_fold_CXX, as.numeric(sub(":.*", "", diff)))
  
  missing<-pdz3_fold_CXX[c(87:89),]
  missing$diff<-c("3:4", "35:36", "20:21")
  missing$insID<-c("RI", "GP", "GE")
  missing$duplicated<-"FALSE"
  
  ##remove
  pdz3_fold_CXX<-pdz3_fold_CXX[-c(87:89),]
  
  ## add the missing variants
  pdz3_fold_CXX<-rbind(pdz3_fold_CXX,
                       missing)
  
  ###re-organise the df based on position
  pdz3_fold_CXX <- arrange(pdz3_fold_CXX, as.numeric(sub(":.*", "", diff)))
  
  ## manually change to the correct positions
  duplicated_pdz3<-pdz3_fold_CXX[c(3,33,58),]
  duplicated_pdz3$diff<-c("3:4", "20:21", "35:36")
  
  pdz3_fold_CXX<-pdz3_fold_CXX[!pdz3_fold_CXX$duplicated == TRUE,]
  
  ## manually change to the correct positions
  duplicated_pdz3$diff<-c("2:3", "21:22", "34:35")
  
  ## add the missing variants
  pdz3_fold_CXX<-rbind(pdz3_fold_CXX,
                       duplicated_pdz3)
  
  ###re-organise the df based on position
  pdz3_fold_CXX <- arrange(pdz3_fold_CXX, as.numeric(sub(":.*", "", diff)))
  
  ## manually change into the correct insertions
  pdz3_fold_CXX[19,]$diff<-"19:20"
  
  ###re-organise the df based on position
  pdz3_fold_CXX <- arrange(pdz3_fold_CXX, as.numeric(sub(":.*", "", diff)))
  
  ## add the missing variants
  missing<-pdz3_fold_CXX[c(4,11,21,23),]
  missing$diff<-c("5:6", "13:14", "24:25", "27:28")
  
  pdz3_fold_CXX<-rbind(pdz3_fold_CXX,
                       missing)
  
  ###re-organise the df based on position
  pdz3_fold_CXX <- arrange(pdz3_fold_CXX, as.numeric(sub(":.*", "", diff)))
  
  ## add mutation type
  pdz3_fold_CXX$type<-"doubleINS_CXX"

  ###re-organise the df based on position
  pdz3_fold_CNN <- arrange(pdz3_fold_CNN, as.numeric(sub(":.*", "", diff)))
  
  ## manually change into the correct insertions and add missing variants
  duplicated_pdz3<-pdz3_fold_CNN[pdz3_fold_CNN$diff %in% c("2:3", "19:20", "34:35"),]
  duplicated_pdz3$diff<-c("3:4", "20:21", "35:36")
  pdz3_fold_CNN<-rbind(pdz3_fold_CNN, 
                       duplicated_pdz3)
  
  ###re-organise the df based on position
  pdz3_fold_CNN <- arrange(pdz3_fold_CNN, as.numeric(sub(":.*", "", diff)))
  
  ## manually change into the correct insertions and add missing variants
  pdz3_fold_CNN<-pdz3_fold_CNN[-c(49:51),]
  missing<-pdz3_fold_insAA[c(15,32),]
  pdz3_fold_CNN<-rbind(pdz3_fold_CNN,
                       missing, fill=T)
  
  ## add mutation type
  pdz3_fold_CNN$type<-"doubleINS_CNN"
  
  
  ### same for the triple ins
  pdz3_fold_CXXX$diff<-as.character(pdz3_fold_CXXX$diff)
  pdz3_fold_CXXX$insID<-substr(pdz3_fold_CXXX$aa_seq, str_first_number(pdz3_fold_CXXX$diff), str_last_number(pdz3_fold_CXXX$diff))
  
  ##now split between CNNN and CXXX by looking for repeated characters in string insID
  pdz3_fold_CXXX$duplicated<- substr(pdz3_fold_CXXX$insID, 1, 2) == substr(pdz3_fold_CXXX$insID, 2, 3)
  
  pdz3_fold_CNNN<-pdz3_fold_CXXX[pdz3_fold_CXXX$duplicated == TRUE,]
  pdz3_fold_CXXX<-pdz3_fold_CXXX[!pdz3_fold_CXXX$duplicated == TRUE,]
  
  ## manually change into the correct insertions and add missing variants
  pdz3_fold_CXXX <- arrange(pdz3_fold_CXXX, as.numeric(sub(":.*", "", diff)))
  
  ## manually change into the correct insertions and add missing variants
  pdz3_fold_CXXX[c(39:45),]$diff<-c("20:22", "24:26", "27:29", "35:37", "13:15", "5:6", "3:4")
  pdz3_fold_CXXX[c(39:45),]$insID<-c("GED", "EGI", "FIS", "GPA", "LGF", "VIH", "RIV")
  pdz3_fold_CXXX <- arrange(pdz3_fold_CXXX, as.numeric(sub(":.*", "", diff)))
  
  missing<-pdz3_fold_CXXX[c(9,19,19,24,24),]
  missing$diff<-c("10:12", "21:23", "22:24", "28:30", "29:31")
  missing$insID<-c("STG", "EDG", "DGE", "ISF", "SFI")
  pdz3_fold_CXXX<-rbind(pdz3_fold_CXXX,
                        missing)
  
  ###re-organise the df based on position
  pdz3_fold_CXXX <- arrange(pdz3_fold_CXXX, as.numeric(sub(":.*", "", diff)))
  pdz3_fold_CNNN <- arrange(pdz3_fold_CNNN, as.numeric(sub(":.*", "", diff)))
  
  ## manually change into the correct insertions and add missing variants
  missing<-pdz3_fold_CNNN[c(2,18,32),]
  missing$diff<-c("3:5", "20:22", "35:37")
  
  pdz3_fold_CNNN<-rbind(pdz3_fold_CNNN,
                        missing)
  
  ###re-organise the df based on position
  pdz3_fold_CNNN <- arrange(pdz3_fold_CNNN, as.numeric(sub(":.*", "", diff)))
  
  
  ## add mutation type
  pdz3_fold_CNNN$type<-"tripleINS_CNNN"
  pdz3_fold_CXXX$type<-"tripleINS_CXXX"
  
  # Return the results
  return(list(pdz3_fold_CX = pdz3_fold_CX,
              pdz3_fold_CXX = pdz3_fold_CXX,
              pdz3_fold_CNN = pdz3_fold_CNN,
              pdz3_fold_insAA = pdz3_fold_insAA,
              pdz3_fold_insCC = pdz3_fold_insCC,
              pdz3_fold_CXXX = pdz3_fold_CXXX,
              pdz3_fold_CNNN = pdz3_fold_CNNN,
              pdz3_fold_singles = pdz3_fold_singles,
              pdz3_fold_allins = pdz3_fold_allins,
              pdz3_fold_Delsub = pdz3_fold_Delsub,
              pdz3_fold_wildtype = pdz3_fold_wildtype,
              pdz3_fold_singleDEL = pdz3_fold_singleDEL,
              pdz3_fold_doubleDEL = pdz3_fold_doubleDEL,
              pdz3_fold_tripleDEL = pdz3_fold_tripleDEL,
              pdz3_fold_synonymous = pdz3_fold_synonymous
  ))
}
  
  