isolate_natural_doms <- function(data_folder){
  
  
  K50_dG_Dataset2_Dataset3 <- read.csv(data_folder)
  
  ##extract the natural domains based on wt pdb name
  K50_dG_Dataset2_Dataset3$pdb<-nchar(K50_dG_Dataset2_Dataset3$WT_name)
  Rocklin_nat_doms<-K50_dG_Dataset2_Dataset3[K50_dG_Dataset2_Dataset3$pdb<9,]
  Rocklin_nat_doms<-Rocklin_nat_doms[!Rocklin_nat_doms$WT_name %like% "con",]
  
  ## add a pdb name column
  Rocklin_nat_doms$pdb_name<-substr(Rocklin_nat_doms$WT_name, 1, 4)
  
  ## isolate all wts
  Rocklin_nat_wt<-Rocklin_nat_doms[Rocklin_nat_doms$mut_type == "wt",]
  Rocklin_nat_wt<-Rocklin_nat_wt[!Rocklin_nat_wt$WT_name %like% "con",]
  
  ##make a data.list of the pdb files
  pdbnames<-unique(Rocklin_nat_doms$WT_name)
  
  ## extract position and mutation
  Rocklin_nat_doms$Pos <- str_first_number(Rocklin_nat_doms$mut_type)
  Rocklin_nat_doms$mut <- NA
  
  Rocklin_nat_doms$mut <- ifelse(substr(Rocklin_nat_doms$mut_type, 1, 1) == "w", "wt",
                                 ifelse(substr(Rocklin_nat_doms$mut_type, 1, 1) == "i", substr(Rocklin_nat_doms$mut_type, 1, 4),
                                        ifelse(substr(Rocklin_nat_doms$mut_type, 1, 1) == "d", "del",
                                               Rocklin_nat_doms$mut_type)))
  
  
  ## split into a list based on pdbname
  data_list_nat_doms <- split(Rocklin_nat_doms, f = Rocklin_nat_doms$WT_name) 
  
  return(list(data_list_nat_doms = data_list_nat_doms,
              pdbnames = pdbnames))
  
}
