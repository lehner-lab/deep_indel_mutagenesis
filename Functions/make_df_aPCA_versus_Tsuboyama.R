make_df_aPCA_versus_Tsuboyama <- function(){
  
  
  single_ins <- scaled_variants_aPCA[scaled_variants_aPCA$type == "singleINS",c("Pos", "domain", "scaled_fitness")]
  colnames(single_ins)[3] <- "scaled_fitness_ins"
  
  single_dels <- scaled_variants_aPCA[scaled_variants_aPCA$type == "singleDEL",c("Pos", "domain", "scaled_fitness")]
  colnames(single_dels)[3] <- "scaled_fitness_dels"
  
  ## align positions to Topolska et al data
  fbp11 <- cor_tsuboyamaVStopolska[cor_tsuboyamaVStopolska$domain == "FBP11-FF1",]
  fbp11 <- fbp11[fbp11$Pos %in% c(9:58),]
  fbp11$Pos <- fbp11$Pos-8
  
  cspacsd <- cor_tsuboyamaVStopolska[cor_tsuboyamaVStopolska$domain == "CSPA-CSD",]
  cspacsd <- cspacsd[cspacsd$Pos %in% c(3:67),]
  cspacsd$Pos <- cspacsd$Pos-2
  
  cspbcsd <- cor_tsuboyamaVStopolska[cor_tsuboyamaVStopolska$domain == "CSPB-CSD",]
  cspbcsd <- cspbcsd[cspbcsd$Pos %in% c(3:65),]
  cspbcsd$Pos <- cspbcsd$Pos-2
  
  bl17ntl9 <- cor_tsuboyamaVStopolska[cor_tsuboyamaVStopolska$domain == "BL17-NTL9",]
  bl17ntl9 <- bl17ntl9[bl17ntl9$Pos %in% c(1:47),]
  
  vil1hp <- cor_tsuboyamaVStopolska[cor_tsuboyamaVStopolska$domain == "VIL1-HP",]
  vil1hp <- vil1hp[vil1hp$Pos %in% c(1:36),]
  
  cor_tsuboyamaVStopolska_realigned <- rbind(fbp11, cspacsd, cspbcsd, bl17ntl9, vil1hp)
  
  cor_tsuboyamaVStopolska_realigned <- merge(cor_tsuboyamaVStopolska_realigned,
                                   single_ins, by = c("Pos", "domain"))
  
  cor_tsuboyamaVStopolska_realigned <- merge(cor_tsuboyamaVStopolska_realigned,
                                   single_dels, by = c("Pos", "domain"))
  

  return(list(cor_tsuboyamaVStopolska_realigned = cor_tsuboyamaVStopolska_realigned))
}

