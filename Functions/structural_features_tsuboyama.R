structural_features_tsuboyama <- function(){
  
  tsuboyama_sec_struc<-bind_rows(list_tsuboyama_nat_doms_strucutre)
  tsuboyama_sec_struc <- tsuboyama_sec_struc[,c("Pos", "pdb_name", "resid", "Structure")]
  tsuboyama_sec_struc <- tsuboyama_sec_struc[!duplicated(tsuboyama_sec_struc),]
  
  list_tsuboyama_sec_struc <- split(tsuboyama_sec_struc, f=tsuboyama_sec_struc$pdb_name)
  
  ## remove the first list as there is only one observation
  list_tsuboyama_sec_struc <- list_tsuboyama_sec_struc[2:length(list_tsuboyama_sec_struc)]
  
  ## first arrange the data by Pos
  for(i in 1:length(list_tsuboyama_sec_struc)) {
    list_tsuboyama_sec_struc[[i]]<-list_tsuboyama_sec_struc[[i]] %>% arrange(Pos)
    
  }
  
  ## translate structure from stride 1-letter code and add termini
  for (i in 1:length(list_tsuboyama_sec_struc)){
    list_tsuboyama_sec_struc[[i]] <- list_tsuboyama_sec_struc[[i]] %>% mutate(Structure=
                                                                                ifelse(Structure == "B", "Bridge",
                                                                                       ifelse(Structure == "C", "Coil",
                                                                                              ifelse(Structure == "T", "Turn", 
                                                                                                     ifelse(Structure == "E", "Strand",
                                                                                                            ifelse(Structure == "G", "310Helix", 
                                                                                                                   ifelse(Structure == "H", "AlphaHelix", NA)))))))
    if (list_tsuboyama_sec_struc[[i]][1,]$Structure %in% c("Coil", "Turn", "Bridge")){
      next_row<-which(!list_tsuboyama_sec_struc[[i]]$Structure %in% c("Coil", "Turn", "Bridge"))
      last_no<-min(next_row)
      list_tsuboyama_sec_struc[[i]][1:(last_no-1),]$Structure<-"ntermini"
    }
  }
  
  for (i in 1:length(list_tsuboyama_sec_struc)){
    if (list_tsuboyama_sec_struc[[i]][nrow(list_tsuboyama_sec_struc[[i]]),]$Structure %in% c("Coil", "Turn", "Bridge")){
      next_row<-which(!list_tsuboyama_sec_struc[[i]]$Structure %in% c("Coil", "Turn", "Bridge"))
      first_no<-max(next_row)
      list_tsuboyama_sec_struc[[i]][(first_no+1):nrow(list_tsuboyama_sec_struc[[i]]),]$Structure <- "ctermini"
    }
  }
  
  ## add simple struc
  for (i in 1:length(list_tsuboyama_sec_struc)){
    list_tsuboyama_sec_struc[[i]] <- list_tsuboyama_sec_struc[[i]] %>% mutate(simple_struc=
                                                                                ifelse(Structure == "Coil" | Structure == "Turn" | Structure == "Bridge", "loop",
                                                                                       ifelse(Structure == "AlphaHelix", "AlphaHelix",
                                                                                              ifelse(Structure == "Strand", "Strand", 
                                                                                                     ifelse(Structure == "310Helix", "310Helix",
                                                                                                            ifelse(Structure == "ntermini", "ntermini",
                                                                                                                   ifelse(Structure == "ctermini", "ctermini", NA)))))))
  }
  
  ## index every secondary structure element in each protein
  # step1: add position for each element
  tsuboyama_sec_struc<-c()
  for(i in 1:length(list_tsuboyama_sec_struc)) {
    temp<-list_tsuboyama_sec_struc[[i]]
    temp$element_pos_simple<-""
    for (i in 1:nrow(temp)){
      if(i==1){
        temp[i,]$element_pos_simple<-"first"
      }else if(temp$simple_struc[i] != temp$simple_struc[i-1]){
        temp[i,]$element_pos_simple<-"first"
      }else if(i==as.numeric(nrow(temp))){
        temp[i,]$element_pos_simple<-"last"
      }else if(temp$simple_struc[i] != temp$simple_struc[i+1]){
        temp[i,]$element_pos_simple<-"last"
      }else if(temp$element_pos_simple[i-1] == "last"){
        temp[i,]$element_pos_simple<-"first"
      }
    }
    for (i in 1:nrow(temp)){
      if(temp[i,]$element_pos_simple == "first"){
        temp[i,]$element_pos_simple<-1
      }else{
        temp[i,]$element_pos_simple<-as.numeric(temp[i-1,]$element_pos_simple)+1
      }
    }
    tsuboyama_sec_struc <- tsuboyama_sec_struc %>%
      rbind(temp)
  }
  
  # step2: add the index
  counter <- 0
  tsuboyama_sec_struc$element_no_simple <- NA
  for (i in 1:nrow(tsuboyama_sec_struc)) {
    if (tsuboyama_sec_struc[i, "element_pos_simple"] == "1") {
      counter <- counter + 1
    }
    tsuboyama_sec_struc[i, "element_no_simple"] <- counter
  }
  
  ## add lenght of each element. 
  tsuboyama_sec_struc$element_no_simple<-as.numeric(tsuboyama_sec_struc$element_no_simple)
  list_tsuboyama_sec_struc <- split(tsuboyama_sec_struc, f = tsuboyama_sec_struc$element_no_simple) 
  
  add_element_lenght<-function(x){
    element_lenght<-as.numeric(nrow(x))
    x<-cbind(x,
             element_lenght)
  }
  
  list_tsuboyama_sec_struc<-lapply(list_tsuboyama_sec_struc,FUN=add_element_lenght)
  
  ## bind rows
  tsuboyama_sec_struc <- bind_rows(list_tsuboyama_sec_struc)
  
  ## add what is before and after each element.
  tsuboyama_sec_struc$structure_before<-""
  tsuboyama_sec_struc$structure_after<-""
  
  ## split by protein
  list_tsuboyama_sec_struc<-split(tsuboyama_sec_struc, f=tsuboyama_sec_struc$pdb_name)
  
  for (i in 1:length(list_tsuboyama_sec_struc)){
    for (j in 1:nrow(list_tsuboyama_sec_struc[[i]])){
      if (j==1){
        list_tsuboyama_sec_struc[[i]][j,]$structure_before<-"start"
      }else if (list_tsuboyama_sec_struc[[i]][j,]$simple_struc != list_tsuboyama_sec_struc[[i]][j-1,]$simple_struc){
        list_tsuboyama_sec_struc[[i]][j,]$structure_before<-list_tsuboyama_sec_struc[[i]][j-1,]$simple_struc
      }else if (list_tsuboyama_sec_struc[[i]][j,]$simple_struc == list_tsuboyama_sec_struc[[i]][j-1,]$simple_struc){
        diff_rows<-which(list_tsuboyama_sec_struc[[i]]$simple_struc != list_tsuboyama_sec_struc[[i]][j,]$simple_struc)
        closest_lower_value<-max(diff_rows[diff_rows<j]) 
        if(is.infinite(closest_lower_value)){
          list_tsuboyama_sec_struc[[i]][j,]$structure_before<-"start"
        }else{
          list_tsuboyama_sec_struc[[i]][j,]$structure_before<-list_tsuboyama_sec_struc[[i]][closest_lower_value,]$simple_struc
        }
      }
    }
  }
  
  
  for (i in 1:length(list_tsuboyama_sec_struc)){
    for (j in 1:nrow(list_tsuboyama_sec_struc[[i]])){
      if (list_tsuboyama_sec_struc[[i]][j,]$simple_struc == list_tsuboyama_sec_struc[[i]][j+1,]$simple_struc & j!=nrow(list_tsuboyama_sec_struc[[i]])){
        diff_rows<-which(list_tsuboyama_sec_struc[[i]]$simple_struc != list_tsuboyama_sec_struc[[i]][j+1,]$simple_struc)
        closest_higher_value<-min(diff_rows[diff_rows>j])
        if (is.infinite(closest_higher_value)){ 
          list_tsuboyama_sec_struc[[i]][j,]$structure_after<-"end"
        }else{
          list_tsuboyama_sec_struc[[i]][j,]$structure_after<-list_tsuboyama_sec_struc[[i]][closest_higher_value,]$simple_struc
        }
      }else if (list_tsuboyama_sec_struc[[i]][j,]$simple_struc != list_tsuboyama_sec_struc[[i]][j+1,]$simple_struc & j!=nrow(list_tsuboyama_sec_struc[[i]])){
        list_tsuboyama_sec_struc[[i]][j,]$structure_after<-list_tsuboyama_sec_struc[[i]][j+1,]$simple_struc
      }else if (j == nrow(list_tsuboyama_sec_struc[[i]]) | is.na(list_tsuboyama_sec_struc[[i]][j+1,]$simple_struc)){
        list_tsuboyama_sec_struc[[i]][j,]$structure_after<-"end"
      }
    }
  }
  
  ## add coil anchor. 
  for (i in 1:length(list_tsuboyama_sec_struc)){
    list_tsuboyama_sec_struc[[i]]$loop_anchor<-""
    for (j in 1:nrow(list_tsuboyama_sec_struc[[i]])){
      list_tsuboyama_sec_struc[[i]][j,]$loop_anchor<-paste(list_tsuboyama_sec_struc[[i]][j,]$structure_before,
                                                           "_",
                                                           list_tsuboyama_sec_struc[[i]][j,]$structure_after,
                                                           sep="")
      
    }
  }
  
  
  ## add secondary structure neighbour (whats neighbouring the secondary structure element).
  secstruct<-c("Strand", "AlphaHelix", "310Helix")
  for (i in 1:length(list_tsuboyama_sec_struc)){
    list_tsuboyama_sec_struc[[i]]$secstruc_neighbour<-""
    
    for (j in 1:nrow(list_tsuboyama_sec_struc[[i]])){
      if (list_tsuboyama_sec_struc[[i]][j,]$simple_struc %in% secstruct){
        current_struc<-list_tsuboyama_sec_struc[[i]][j,]$simple_struc
        current_no<-list_tsuboyama_sec_struc[[i]][j,]$element_no_simple
        
        neighbour<-which(list_tsuboyama_sec_struc[[i]]$simple_struc != "loop"
                         & list_tsuboyama_sec_struc[[i]]$element_no_simple != current_no)
        
        before<-max(neighbour[neighbour<j])
        
        after<-min(neighbour[neighbour>j])
        
        if (!is.infinite(before) & !is.infinite(after)){
          list_tsuboyama_sec_struc[[i]][j,]$secstruc_neighbour<-paste(list_tsuboyama_sec_struc[[i]][before,]$simple_struc,
                                                                      "_",
                                                                      list_tsuboyama_sec_struc[[i]][after,]$simple_struc,
                                                                      sep="")
        }else if (is.infinite(before) & !is.infinite(after)){
          list_tsuboyama_sec_struc[[i]][j,]$secstruc_neighbour<-paste("start",
                                                                      "_",
                                                                      list_tsuboyama_sec_struc[[i]][after,]$simple_struc,
                                                                      sep="")
          
        }else if (!is.infinite(before) & is.infinite(after)){
          list_tsuboyama_sec_struc[[i]][j,]$secstruc_neighbour<-paste(list_tsuboyama_sec_struc[[i]][before,]$simple_struc,
                                                                      "_",
                                                                      "end",
                                                                      sep="")
        }
      }
    }
  }
  
  
  tsuboyama_sec_struc<-bind_rows(list_tsuboyama_sec_struc)
  
  ## add the rsasa: tsuboyama_rsasa
  tsuboyama_sec_struc <- merge(tsuboyama_sec_struc,
                               tsuboyama_rsasa,
                               by = c("Pos", "pdb_name"),
                               all.x = T)
  
  ## realign the secondary structure using rSASA
  # the idea is to realign positions for each secondary structure element seperatly counting from middle to both terminis
  # step 1: find median of secondary element. 
  # step 2: set position 0 to the position that is closest to the core of the proteins (<rSASA), but also only +/-1 aa away from median. 
  
  # split based on element
  list_tsuboyama_sec_struc <- split(tsuboyama_sec_struc, f = tsuboyama_sec_struc$element_no_simple) 
  
  ##sort based on the element_pos_simple
  for(i in 1:length(list_tsuboyama_sec_struc)) {
    list_tsuboyama_sec_struc[[i]]$element_pos_simple<-as.numeric(list_tsuboyama_sec_struc[[i]]$element_pos_simple)
    list_tsuboyama_sec_struc[[i]]<-list_tsuboyama_sec_struc[[i]] %>% arrange(element_pos_simple)
    
  }
  
  # this script carries out step 1 and 2 for each secondary structure element. 
  for (i in 1:length(list_tsuboyama_sec_struc)) {
    list_tsuboyama_sec_struc[[i]]$element_pos_simple<-as.numeric(list_tsuboyama_sec_struc[[i]]$element_pos_simple)
    median<-median(list_tsuboyama_sec_struc[[i]]$element_pos_simple)
    if (nrow(list_tsuboyama_sec_struc[[i]])>2){
      if (unique(list_tsuboyama_sec_struc[[i]]$simple_struc) == "AlphaHelix"){
        if (median %% 2 %in% c("0","1")){
          if (list_tsuboyama_sec_struc[[i]][median,]$rSASA < list_tsuboyama_sec_struc[[i]][median+1,]$rSASA & list_tsuboyama_sec_struc[[i]][median,]$rSASA < list_tsuboyama_sec_struc[[i]][median-1,]$rSASA){
            median<-median
          }else if(list_tsuboyama_sec_struc[[i]][median,]$rSASA > list_tsuboyama_sec_struc[[i]][median+1,]$rSASA & list_tsuboyama_sec_struc[[i]][median,]$rSASA < list_tsuboyama_sec_struc[[i]][median-1,]$rSASA){
            median<-median+1
          }else if(list_tsuboyama_sec_struc[[i]][median,]$rSASA < list_tsuboyama_sec_struc[[i]][median+1,]$rSASA & list_tsuboyama_sec_struc[[i]][median,]$rSASA > list_tsuboyama_sec_struc[[i]][median-1,]$rSASA){
            median<-median-1
          }else if(list_tsuboyama_sec_struc[[i]][median,]$rSASA > list_tsuboyama_sec_struc[[i]][median+1,]$rSASA & list_tsuboyama_sec_struc[[i]][median,]$rSASA > list_tsuboyama_sec_struc[[i]][median-1,]$rSASA & nrow(list_tsuboyama_sec_struc[[i]])>4){
            median<-median+2
          }
        }else if(median %% 2 %in% c("0.5","1.5")){
          median<-median-0.5
          if (list_tsuboyama_sec_struc[[i]][median,]$rSASA < list_tsuboyama_sec_struc[[i]][median+1,]$rSASA & list_tsuboyama_sec_struc[[i]][median,]$rSASA < list_tsuboyama_sec_struc[[i]][median-1,]$rSASA){
            median<-median
          }else if(list_tsuboyama_sec_struc[[i]][median,]$rSASA > list_tsuboyama_sec_struc[[i]][median+1,]$rSASA & list_tsuboyama_sec_struc[[i]][median,]$rSASA < list_tsuboyama_sec_struc[[i]][median-1,]$rSASA){
            median<-median+1
          }else if(list_tsuboyama_sec_struc[[i]][median,]$rSASA < list_tsuboyama_sec_struc[[i]][median+1,]$rSASA & list_tsuboyama_sec_struc[[i]][median,]$rSASA > list_tsuboyama_sec_struc[[i]][median-1,]$rSASA){
            median<-median-1
          }else if(list_tsuboyama_sec_struc[[i]][median,]$rSASA > list_tsuboyama_sec_struc[[i]][median+1,]$rSASA & list_tsuboyama_sec_struc[[i]][median,]$rSASA > list_tsuboyama_sec_struc[[i]][median-1,]$rSASA & nrow(list_tsuboyama_sec_struc[[i]])>4){
            median<-median+2
          }
        }
      }else if (unique(list_tsuboyama_sec_struc[[i]]$simple_struc) == "310Helix"){
        if (nrow(list_tsuboyama_sec_struc[[i]]) %in% c(3,4)){
          median<-which.min(list_tsuboyama_sec_struc[[i]]$rSASA)
        }else if (nrow(list_tsuboyama_sec_struc[[i]])>4){
          median<-which.min(list_tsuboyama_sec_struc[[i]]$rSASA)
        }
      }else if (unique(list_tsuboyama_sec_struc[[i]]$simple_struc) == "Strand"){
        if (median %% 2 %in% c("0","1")){
          if (list_tsuboyama_sec_struc[[i]][median,]$rSASA > list_tsuboyama_sec_struc[[i]][median+1,]$rSASA & list_tsuboyama_sec_struc[[i]][median,]$rSASA > list_tsuboyama_sec_struc[[i]][median-1,]$rSASA){
            median<-median+1
          }else if(list_tsuboyama_sec_struc[[i]][median,]$rSASA < list_tsuboyama_sec_struc[[i]][median+1,]$rSASA & list_tsuboyama_sec_struc[[i]][median,]$rSASA < list_tsuboyama_sec_struc[[i]][median-1,]$rSASA){
            median<-median
          }
        }else if(median %% 2 %in% c("0.5","1.5")){
          median<-median-0.5
          if (list_tsuboyama_sec_struc[[i]][median,]$rSASA > list_tsuboyama_sec_struc[[i]][median+1,]$rSASA & list_tsuboyama_sec_struc[[i]][median,]$rSASA > list_tsuboyama_sec_struc[[i]][median-1,]$rSASA){
            median<-median+1
          }else if(list_tsuboyama_sec_struc[[i]][median,]$rSASA < list_tsuboyama_sec_struc[[i]][median+1,]$rSASA & list_tsuboyama_sec_struc[[i]][median,]$rSASA < list_tsuboyama_sec_struc[[i]][median-1,]$rSASA){
            median<-median
          }
        }
      }else if (unique(list_tsuboyama_sec_struc[[i]]$simple_struc) == "loop"){
        if (median %% 2 %in% c("0","1")){
          median<-median
        }else if(median %% 2 %in% c("0.5","1.5")){
          median<-median-0.5
        }
      }else if (unique(list_tsuboyama_sec_struc[[i]]$simple_struc) == "ntermini"){
        median<-nrow(list_tsuboyama_sec_struc[[i]])
      }else if (unique(list_tsuboyama_sec_struc[[i]]$simple_struc) == "ctermini"){
        median<-1
      }
      if (unique(list_tsuboyama_sec_struc[[i]]$simple_struc) %in% c("AlphaHelix", "310Helix", "Strand", "loop")){
        list_tsuboyama_sec_struc[[i]]$align_to_center<-""
        list_tsuboyama_sec_struc[[i]][median,]$align_to_center<-0
        if (list_tsuboyama_sec_struc[[i]][nrow(list_tsuboyama_sec_struc[[i]]),]$align_to_center == "0"){
          list_tsuboyama_sec_struc[[i]][(nrow(list_tsuboyama_sec_struc[[i]])-1):1,]$align_to_center<- -1:-(median-1)
        }else if (list_tsuboyama_sec_struc[[i]][1,]$align_to_center == "0"){
          list_tsuboyama_sec_struc[[i]][(median+1):nrow(list_tsuboyama_sec_struc[[i]]),]$align_to_center<-1:(nrow(list_tsuboyama_sec_struc[[i]])-1)
        }else{
          list_tsuboyama_sec_struc[[i]][(median+1):nrow(list_tsuboyama_sec_struc[[i]]),]$align_to_center<-c(1:nrow(list_tsuboyama_sec_struc[[i]][(median+1):nrow(list_tsuboyama_sec_struc[[i]]),]))
          list_tsuboyama_sec_struc[[i]][(median-1):1,]$align_to_center<-c(-1:-(nrow(list_tsuboyama_sec_struc[[i]][(median-1):1,])))
        }
      }else if (unique(list_tsuboyama_sec_struc[[i]]$simple_struc) == "ntermini"){
        list_tsuboyama_sec_struc[[i]]$align_to_center<-""
        list_tsuboyama_sec_struc[[i]][c(1:(nrow(list_tsuboyama_sec_struc[[i]]))),]$align_to_center<-c(-(nrow(list_tsuboyama_sec_struc[[i]])):-1)
      }else if (unique(list_tsuboyama_sec_struc[[i]]$simple_struc) == "ctermini"){
        list_tsuboyama_sec_struc[[i]]$align_to_center<-""
        list_tsuboyama_sec_struc[[i]][c(1:(nrow(list_tsuboyama_sec_struc[[i]]))),]$align_to_center<-c(1:(nrow(list_tsuboyama_sec_struc[[i]])))
      }
    }else if (nrow(list_tsuboyama_sec_struc[[i]]) == 1 & !unique(list_tsuboyama_sec_struc[[i]]$simple_struc) %in% c("ntermini", "ctermini")){
      list_tsuboyama_sec_struc[[i]]$align_to_center<-""
      list_tsuboyama_sec_struc[[i]]$align_to_center<-0
    }else if (nrow(list_tsuboyama_sec_struc[[i]]) == 2 & !unique(list_tsuboyama_sec_struc[[i]]$simple_struc) %in% c("ntermini", "ctermini")){
      if (list_tsuboyama_sec_struc[[i]][1,]$rSASA < list_tsuboyama_sec_struc[[i]][2,]$rSASA){
        list_tsuboyama_sec_struc[[i]]$align_to_center<-""
        list_tsuboyama_sec_struc[[i]]$align_to_center<-c(0,1)
      }else{
        list_tsuboyama_sec_struc[[i]]$align_to_center<-""
        list_tsuboyama_sec_struc[[i]]$align_to_center<-c(-1,0)
      }
    }else if (nrow(list_tsuboyama_sec_struc[[i]]) == 1 & unique(list_tsuboyama_sec_struc[[i]]$simple_struc) %in% c("ntermini")){
      list_tsuboyama_sec_struc[[i]]$align_to_center<-""
      list_tsuboyama_sec_struc[[i]]$align_to_center<-c(-1)
    }else if (nrow(list_tsuboyama_sec_struc[[i]]) == 1 & unique(list_tsuboyama_sec_struc[[i]]$simple_struc) %in% c("ctermini")){
      list_tsuboyama_sec_struc[[i]]$align_to_center<-""
      list_tsuboyama_sec_struc[[i]]$align_to_center<-c(1)
    }else if (nrow(list_tsuboyama_sec_struc[[i]]) == 2 & unique(list_tsuboyama_sec_struc[[i]]$simple_struc) %in% c("ntermini")){
      list_tsuboyama_sec_struc[[i]]$align_to_center<-""
      list_tsuboyama_sec_struc[[i]]$align_to_center<-c(-1,-2)
    }else if (nrow(list_tsuboyama_sec_struc[[i]]) == 2 & unique(list_tsuboyama_sec_struc[[i]]$simple_struc) %in% c("ctermini")){
      list_tsuboyama_sec_struc[[i]]$align_to_center<-""
      list_tsuboyama_sec_struc[[i]]$align_to_center<-c(1,2)
    }
  }
  
  ## convert the new position metric "align_to_centre" to numeric
  for (i in 1:length(list_tsuboyama_sec_struc)){
    list_tsuboyama_sec_struc[[i]]$align_to_center<-as.numeric(list_tsuboyama_sec_struc[[i]]$align_to_center)
  }
  
  tsuboyama_sec_struc<-bind_rows(list_tsuboyama_sec_struc)
  
  ## add a column for secondary structure elements that are directly followed by another structure element. 
  list_tsuboyama_sec_struc <- split(tsuboyama_sec_struc, f=tsuboyama_sec_struc$pdb_name)
  secstruct<-c("Strand", "AlphaHelix", "310Helix")
  
  for (i in 1:length(list_tsuboyama_sec_struc)){
    list_tsuboyama_sec_struc[[i]]$neighbour_n<-""
    list_tsuboyama_sec_struc[[i]]$neighbour_c<-""
    for (j in 1:nrow(list_tsuboyama_sec_struc[[i]])){
      if (list_tsuboyama_sec_struc[[i]][j,]$simple_struc %in% secstruct){
        list_tsuboyama_sec_struc[[i]][j,]$neighbour_n<-gsub("_.*", "", list_tsuboyama_sec_struc[[i]][j,]$secstruc_neighbour)
        list_tsuboyama_sec_struc[[i]][j,]$neighbour_c<-gsub(".*_", "", list_tsuboyama_sec_struc[[i]][j,]$secstruc_neighbour)
        
      }
    }
  }
  
  ##  remove the loop_anchor if the simple_struc is not loop
  for (i in 1:length(list_tsuboyama_sec_struc)){
    for (j in 1:nrow(list_tsuboyama_sec_struc[[i]])){
      if (list_tsuboyama_sec_struc[[i]][j,]$simple_struc %in% secstruct){
        list_tsuboyama_sec_struc[[i]][j,]$loop_anchor<-""
      }
    }
  }
  
  ## collpase final df
  tsuboyama_sec_struc<-bind_rows(list_tsuboyama_sec_struc)
  
  return(list(tsuboyama_sec_struc = tsuboyama_sec_struc))
  
}












