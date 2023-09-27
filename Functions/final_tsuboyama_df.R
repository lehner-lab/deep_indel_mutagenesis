final_tsuboyama_df <- function(){
  
  ## collapse list: list_tsuboyama_nat_doms_strucutre
  tsuboyama_nat_doms_structure<-bind_rows(list_tsuboyama_nat_doms_strucutre)
  
  ## remove unneccesery columns
  tsuboyama_nat_doms_structure<-tsuboyama_nat_doms_structure[,c("Pos", "aa_seq", "WT_name", 
                                                                "ddG_ML", "pdb_name", "mut", 
                                                                "resid", "Structure")]
  
  ## split df into insA, insG, dels and single subs.
  tsuboyama_nat_doms_structure_insA<-tsuboyama_nat_doms_structure[tsuboyama_nat_doms_structure$mut == "insA",]
  tsuboyama_nat_doms_structure_insG<-tsuboyama_nat_doms_structure[tsuboyama_nat_doms_structure$mut == "insG",]
  
  tsuboyama_nat_doms_structure_dels<-tsuboyama_nat_doms_structure[tsuboyama_nat_doms_structure$mut %like% "del",]
  
  tsuboyama_nat_doms_structure_subs<-tsuboyama_nat_doms_structure[!tsuboyama_nat_doms_structure$mut %like% "del",]
  tsuboyama_nat_doms_structure_subs<-tsuboyama_nat_doms_structure_subs[!tsuboyama_nat_doms_structure_subs$mut %like% "ins",]
  tsuboyama_nat_doms_structure_subs<-tsuboyama_nat_doms_structure_subs[!tsuboyama_nat_doms_structure_subs$mut %like% ":",]
  
  
  ## calculate mean sub ddG/position
  tsuboyama_nat_doms_structure_subs$ddG_ML <- as.numeric(tsuboyama_nat_doms_structure_subs$ddG_ML)
  tsuboyama_nat_doms_structure_subs_avg <- aggregate(tsuboyama_nat_doms_structure_subs$ddG_ML, list(tsuboyama_nat_doms_structure_subs$Pos, tsuboyama_nat_doms_structure_subs$pdb_name), FUN=mean)
  colnames(tsuboyama_nat_doms_structure_subs_avg)[1:3]<-c("Pos", "pdb_name", "ddG_ML_avg")

  ### now create a merge data set for all mutation types. 
  tsuboyama_nat_doms_indels<-merge(tsuboyama_nat_doms_structure_dels[,c("Pos", "ddG_ML", "pdb_name", "resid", "Structure")],
                                   tsuboyama_nat_doms_structure_insA[,c("Pos", "ddG_ML", "pdb_name")],
                                   by=c("pdb_name","Pos"))
  
  colnames(tsuboyama_nat_doms_indels)[3] <- c("ddG_ML_dels")
  colnames(tsuboyama_nat_doms_indels)[6] <- c("ddG_ML_ins")
  
  tsuboyama_nat_doms_all<-merge(tsuboyama_nat_doms_indels,
                                tsuboyama_nat_doms_structure_subs_avg,
                                by=c("Pos", "pdb_name"))
  
  colnames(tsuboyama_nat_doms_all)[7]<-c("ddG_ML_subs")
  
  ## remove ddG for ins, dels or subs that are NA. 
  tsuboyama_nat_doms_all$ddG_ML_dels <- as.numeric(tsuboyama_nat_doms_all$ddG_ML_dels)
  tsuboyama_nat_doms_all$ddG_ML_ins <- as.numeric(tsuboyama_nat_doms_all$ddG_ML_ins)
  
  tsuboyama_nat_doms_all<-tsuboyama_nat_doms_all[-which(is.na(tsuboyama_nat_doms_all$ddG_ML_dels)),]
  tsuboyama_nat_doms_all<-tsuboyama_nat_doms_all[-which(is.na(tsuboyama_nat_doms_all$ddG_ML_ins)),]
  tsuboyama_nat_doms_all<-tsuboyama_nat_doms_all[-which(is.na(tsuboyama_nat_doms_all$ddG_ML_subs)),]
  
  ## reverse ddG sign
  tsuboyama_nat_doms_all$ddG_ML_dels <- tsuboyama_nat_doms_all$ddG_ML_dels* -1
  tsuboyama_nat_doms_all$ddG_ML_ins <- tsuboyama_nat_doms_all$ddG_ML_ins* -1
  tsuboyama_nat_doms_all$ddG_ML_subs <- tsuboyama_nat_doms_all$ddG_ML_subs* -1
  tsuboyama_nat_doms_structure_subs$ddG_ML <- tsuboyama_nat_doms_structure_subs$ddG_ML* -1
  
  
  return(list(tsuboyama_nat_doms_all = tsuboyama_nat_doms_all,
              tsuboyama_nat_doms_structure_subs = tsuboyama_nat_doms_structure_subs))
  
}
