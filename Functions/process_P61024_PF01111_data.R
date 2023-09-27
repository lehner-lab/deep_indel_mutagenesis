
### this script defines the function to call the different mutational variants for P61024_PF01111
process_P61024_PF01111_data <- function(input_file_path) {
  # Load the RData file
  load(input_file_path)
  
  # Create a copy of the data for the domain
  P61024_PF01111 <- all_variants
  
  # Identify positions where sequences differ from wt_seq at the aa level
  P61024_PF01111$diff <- mapply(
    function(x, y) grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "IYYSDKYDDEEFEYRHVMLPKDIAKLVPKTHLMSESEWRNLGVQQSQGWVHYMIHEPEPHILLFRRP",
    P61024_PF01111$aa_seq
  )
  
  # Add length of the sequence to the data frame
  P61024_PF01111$lenght <- str_length(P61024_PF01111$aa_seq)
  
  # Count the number of characters where sequences differ from wt_seq
  P61024_PF01111$characters <- nchar(P61024_PF01111$diff)
  
  # Isolate the correct sequences
  P61024_PF01111_variants <- P61024_PF01111[P61024_PF01111$characters < 14, ]
  
  # Find correct synonymous mutants for this wt
  P61024_PF01111_syns <- synonymous
  
  # Identify positions where sequences differ from wt_seq at the nt level
  P61024_PF01111_syns$diff <- mapply(
    function(x, y) grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "atttattattctgataaatatgatgatgaagaatttgaatatagacatgttatgttgccaaaagatattgcaaaattggttccaaaaactcatttgatgtctgaatctgaatggagaaatttgggtgttcaacaatctcaaggttgggttcattatatgattcatgaaccagaaccacatattttgttgtttagaagacca",
    P61024_PF01111_syns$nt_seq
  )
  
  # Count the number of characters where sequences differ from wt_seq
  P61024_PF01111_syns$characters <- nchar(P61024_PF01111_syns$diff)
  
  # Filter the data frame for synonymous mutants
  P61024_PF01111_syns <- P61024_PF01111_syns[P61024_PF01111_syns$characters < 12, ]
  
  
  
  # Return the results
  return(list(P61024_PF01111_syns = P61024_PF01111_syns, 
              P61024_PF01111_variants = P61024_PF01111_variants))
}

