

find_subID_delsubs_grb2 <- function(){

  ## isolate delsubs
  grb2_delsubs <- scaled_variants_bPCA[scaled_variants_bPCA$domain == "GRB2-SH3" &
                                       scaled_variants_bPCA$mut_type == "delSub",]
  grb2_delsubs$subID <- substr(grb2_delsubs$aa_seq, grb2_delsubs$Pos, grb2_delsubs$Pos)
  
  grb2_delsubs_1 <- grb2_delsubs[grb2_delsubs$type == "delSub_1",]
  grb2_delsubs_2 <- grb2_delsubs[grb2_delsubs$type == "delSub_2",]
  
  ## make a df to fill out the missing values. 
  missing_rows<-data.table(Pos=c(1:52),
                         subID="-")
  
  rows_to_keep<-which(!missing_rows$Pos %in% grb2_delsubs_1$Pos)
  missing_rows <- missing_rows[rows_to_keep,]
  
  grb2_delsubs_1 <- rbind(grb2_delsubs_1,
                        missing_rows, fill =T)
  
  missing_rows<-data.table(Pos=c(1:52),
                         subID="-")
  
  rows_to_keep<-which(!missing_rows$Pos %in% grb2_delsubs_2$Pos)
  missing_rows <- missing_rows[rows_to_keep,]
  
  grb2_delsubs_2 <- rbind(grb2_delsubs_2,
                        missing_rows, fill =T)

  return(list(grb2_delsubs_1 = grb2_delsubs_1,
              grb2_delsubs_2 = grb2_delsubs_2))

}
