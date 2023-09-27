run_stride <- function(data_folder){
  
  ## read pdb files
  setwd(data_folder)
  
  pdb_files <- paste(pdbnames, 
                     sep="")
  
  read_pdb_files <- lapply(pdb_files, read.pdb, verbose=F)
  
  ## run stride
  stride_files<-lapply(read_pdb_files, stride, resno=TRUE)
  
  ## extract pdb info info into a df
  extract_pdb_info<-function(x){
    x<-data.frame(x$atom)
    x<-x[,c(5:7)]
    x<- x %>% distinct(resno, chain, .keep_all=T)
    x$resid<-a(paste(substr(x$resid,1,1),
                     tolower(substr(x$resid,2,3)),
                     sep=""))
    
    x<-x[!is.na(x$resid),]
    return(x)
  }
  
  strucutre_info_tsuboyama<-lapply(read_pdb_files,FUN=extract_pdb_info)
  
  ## extract stride info info into a df
  extract_stride_info<-function(y){
    
    y<-data.frame(Structure=y$sse, ACC=y$acc)
    return(y)
  }
  
  strucutre_info_tsuboyama_2<-lapply(stride_files,FUN=extract_stride_info)
  
  ### now merge pdb and stride info. 
  strucutre_tsuboyama_list<- NULL
  for (i in 1:length(strucutre_info_tsuboyama)) {
    strucutre_tsuboyama_list[[i]]<-bind_cols(strucutre_info_tsuboyama[i], strucutre_info_tsuboyama_2[i])
  }
  
  ## rename after pdbname
  names(strucutre_tsuboyama_list)<-paste(toupper(substr(pdb_files, 1, 4)),
                                       ".pdb",
                                       sep="")
  
  ## rename resno to pos
  for (i in 1:length(strucutre_tsuboyama_list)) {
    setnames(strucutre_tsuboyama_list[[i]], old = "resno", 
             new = "Pos")
  }
  
  return(strucutre_tsuboyama_list = strucutre_tsuboyama_list)
  
}


