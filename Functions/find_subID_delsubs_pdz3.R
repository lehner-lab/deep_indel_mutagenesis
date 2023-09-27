
find_subID_delsubs_pdz3 <- function(){
  
  ## isolate delsubs
  pdz3_delsubs <- scaled_variants_bPCA[scaled_variants_bPCA$domain == "PSD95-PDZ3" &
                                         scaled_variants_bPCA$mut_type == "delSub",]
  pdz3_delsubs$subID <- substr(pdz3_delsubs$aa_seq, pdz3_delsubs$Pos, pdz3_delsubs$Pos)
  
  pdz3_delsubs_1 <- pdz3_delsubs[pdz3_delsubs$type == "delSub_1",]
  pdz3_delsubs_2 <- pdz3_delsubs[pdz3_delsubs$type == "delSub_2",]
  
  ## make a df to fill out the missing values. 
  missing_rows<-data.table(Pos=c(1:52),
                           subID="-")
  
  rows_to_keep<-which(!missing_rows$Pos %in% pdz3_delsubs_1$Pos)
  missing_rows <- missing_rows[rows_to_keep,]
  
  pdz3_delsubs_1 <- rbind(pdz3_delsubs_1,
                          missing_rows, fill =T)
  
  missing_rows<-data.table(Pos=c(1:52),
                           subID="-")
  
  rows_to_keep<-which(!missing_rows$Pos %in% pdz3_delsubs_2$Pos)
  missing_rows <- missing_rows[rows_to_keep,]
  
  pdz3_delsubs_2 <- rbind(pdz3_delsubs_2,
                          missing_rows, fill =T)
  
  return(list(pdz3_delsubs_1 = pdz3_delsubs_1,
              pdz3_delsubs_2 = pdz3_delsubs_2))
  
}
