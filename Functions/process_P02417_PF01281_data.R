
### this script defines the function to call the different mutational variants for P02417_PF01281
process_P02417_PF01281_data <- function(input_file_path) {
  # Load the RData file
  load(input_file_path)
  
  # Create a copy of the data for the domain
  P02417_PF01281 <- all_variants
  
  # Identify positions where sequences differ from wt_seq at the aa level
  P02417_PF01281$diff <- mapply(
    function(x, y) grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "MKVIFLKDVKGKGKKGEIKNVADGYANNFLFKQGLAIEATPANLKAL",
    P02417_PF01281$aa_seq
  )
  
  # Add length of the sequence to the data frame
  P02417_PF01281$lenght <- str_length(P02417_PF01281$aa_seq)
  
  # Count the number of characters where sequences differ from wt_seq
  P02417_PF01281$characters <- nchar(P02417_PF01281$diff)
  
  # Isolate the correct sequences
  P02417_PF01281_variants <- P02417_PF01281[P02417_PF01281$characters < 14, ]
  
  # Find correct synonymous mutants for this wt
  P02417_PF01281_syns <- synonymous
  
  # Identify positions where sequences differ from wt_seq at the nt level
  P02417_PF01281_syns$diff <- mapply(
    function(x, y) grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "atgaaagttatttttttgaaagatgttaaaggtaaaggtaaaaaaggtgaaattaaaaatgttgcagatggttatgcaaataattttttgtttaaacaaggtttggcaattgaagcaactccagcaaatttgaaagcattg",
    P02417_PF01281_syns$nt_seq
  )
  
  # Count the number of characters where sequences differ from wt_seq
  P02417_PF01281_syns$characters <- nchar(P02417_PF01281_syns$diff)
  
  # Filter the data frame for synonymous mutants
  P02417_PF01281_syns <- P02417_PF01281_syns[P02417_PF01281_syns$characters < 12, ]
  
  
  
  # Return the results
  return(list(P02417_PF01281_syns = P02417_PF01281_syns, 
              P02417_PF01281_variants = P02417_PF01281_variants))
}

