load_deletion_prediction_data <- function(output_folder){
  
  
  file_names <- list.files(output_folder, pattern = "cross_val_\\w+\\.txt", full.names = TRUE)
  
  # Import files into a list
  file_list <- lapply(file_names, read.table)
  
  # Import files into a list with column names
  file_list <- lapply(file_names, function(file) {
    data <- read.table(file, header = TRUE)  # Set header = TRUE to use the first row as column names
    return(data)
  })
  
  ### add names to the list
  list_of_names<-c()
  for (i in 1:length(file_list)){
    domain<-unique(file_list[[i]]$domain)
    
    list_of_names<-append(list_of_names,domain)
  }

  names(file_list) <- list_of_names
  
  return(file_list = file_list)
  
}
