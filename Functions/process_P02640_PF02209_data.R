
### this script defines the function to call the different mutational variants for P02640_PF02209
process_P02640_PF02209_data <- function(input_file_path) {
  # Load the RData file
  load(input_file_path)
  
  # Create a copy of the data for the domain
  P02640_PF02209 <- all_variants
  
  # Identify positions where sequences differ from wt_seq at the aa level
  P02640_PF02209$diff <- mapply(
    function(x, y) grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "HLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF",
    P02640_PF02209$aa_seq
  )
  
  # Add length of the sequence to the data frame
  P02640_PF02209$lenght <- str_length(P02640_PF02209$aa_seq)
  
  # Count the number of characters where sequences differ from wt_seq
  P02640_PF02209$characters <- nchar(P02640_PF02209$diff)
  
  # Isolate the correct sequences
  P02640_PF02209_variants <- P02640_PF02209[P02640_PF02209$characters < 14, ]
  
  # Find correct synonymous mutants for this wt
  P02640_PF02209_syns <- synonymous
  
  # Identify positions where sequences differ from wt_seq at the nt level
  P02640_PF02209_syns$diff <- mapply(
    function(x, y) grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "catttgtctgatgaagattttaaagcagtttttggtatgactagatctgcatttgcaaatttgccattgtggaaacaacaaaatttgaaaaaagaaaaaggtttgttt",
    P02640_PF02209_syns$nt_seq
  )
  
  # Count the number of characters where sequences differ from wt_seq
  P02640_PF02209_syns$characters <- nchar(P02640_PF02209_syns$diff)
  
  # Filter the data frame for synonymous mutants
  P02640_PF02209_syns <- P02640_PF02209_syns[P02640_PF02209_syns$characters < 12, ]
  
  
  
  # Return the results
  return(list(P02640_PF02209_syns = P02640_PF02209_syns, 
              P02640_PF02209_variants = P02640_PF02209_variants))
}

