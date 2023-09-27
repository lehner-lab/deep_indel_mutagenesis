make_pdb_cspacsd_figure1 <- function(x, y, z, q){
  
  ## read the pdb
  pdb<-read.pdb(x)
  
  # make a df of mean insertion/position
  mean_ins <- data.table(aggregate(scaled_variants_aPCA[scaled_variants_aPCA$domain == y &
                                                          scaled_variants_aPCA$mut_type == "insertions",]$scaled_fitness, list(scaled_variants_aPCA[scaled_variants_aPCA$domain == y &
                                                                                                                                                      scaled_variants_aPCA$mut_type == "insertions",]$Pos), FUN = mean))
  colnames(mean_ins)<-c("position","mean_fitness")
  
  ## add the number to align positions to the positions defined in the pdb file
  mean_ins$position<-mean_ins$position+z
  
  
  #set b-factor to mutation data
  unique(pdb$atom$resno[pdb$atom$chain=="A"])
  
  for(i in 5:69){
    pdb$atom$b[pdb$atom$resno==i & pdb$atom$chain=="A"]<-mean_ins[position==i,]$mean_fitness
  }
  
  #set b-factor to zero outside mutated region
  for(i in c(1:4,67)){
    pdb$atom$b[pdb$atom$resno==i & pdb$atom$chain=="A"]<-1
  }
  
  write.pdb(pdb,file=paste(q,
                           "/",
                           y,
                           "_insertions_mean.pdb",
                           sep=""))
  
  ## now for deletions
  pdb<-read.pdb(x)
  
  # make a df of mean insertion/position
  mean_dels <- data.table(aggregate(scaled_variants_aPCA[scaled_variants_aPCA$domain == y &
                                                           scaled_variants_aPCA$mut_type == "deletions",]$scaled_fitness, list(scaled_variants_aPCA[scaled_variants_aPCA$domain == y &
                                                                                                                                                      scaled_variants_aPCA$mut_type == "deletions",]$Pos), FUN = mean))
  colnames(mean_dels)<-c("position","mean_fitness")
  
  ## add the number to align positions to the positions defined in the pdb file
  mean_dels$position<-mean_dels$position+z
  
  #set b-factor to mutation data
  unique(pdb$atom$resno[pdb$atom$chain=="A"])
  for(i in 5:69){
    pdb$atom$b[pdb$atom$resno==i & pdb$atom$chain=="A"]<-mean_dels[position==i,]$mean_fitness
  }
  
  #set b-factor to zero outside mutated region
  for(i in c(1:4,67)){
    pdb$atom$b[pdb$atom$resno==i & pdb$atom$chain=="A"]<-1
  }
  
  write.pdb(pdb,file=paste(q,
                           "/",
                           y,
                           "_deletions_mean.pdb",
                           sep=""))
  
  
}
