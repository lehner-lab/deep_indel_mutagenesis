
### this script defines the function to call the different mutational variants for P01053_PF00280
process_P01053_PF00280_data <- function(input_file_path) {
  # Load the RData file
  load(input_file_path)
  
  # Create a copy of the data for the domain
  P01053_PF00280 <- all_variants
  
  # Identify positions where sequences differ from wt_seq at the aa level
  P01053_PF00280$diff <- mapply(
    function(x, y) grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "KTEWPELVGKSVEEAKKVILQDKPEAQIIVLPVGTIVTMEYRIDRVRLFVDKLDNIAQVPRVG",
    P01053_PF00280$aa_seq
  )
  
  # Add length of the sequence to the data frame
  P01053_PF00280$lenght <- str_length(P01053_PF00280$aa_seq)
  
  # Count the number of characters where sequences differ from wt_seq
  P01053_PF00280$characters <- nchar(P01053_PF00280$diff)
  
  # Isolate the correct sequences
  P01053_PF00280_variants <- P01053_PF00280[P01053_PF00280$characters < 14, ]
  
  # Find correct synonymous mutants for this wt
  P01053_PF00280_syns <- synonymous
  
  # Identify positions where sequences differ from wt_seq at the nt level
  P01053_PF00280_syns$diff <- mapply(
    function(x, y) grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "aaaactgaatggccagaattggttggtaaatctgttgaagaagcaaaaaaagttattttgcaagataaaccagaagcacaaattattgttttgccagttggtactattgttactatggaatatagaattgatagagttagattgtttgttgataaattggataatattgcacaagttccaagagttggt",
    P01053_PF00280_syns$nt_seq
  )
  
  # Count the number of characters where sequences differ from wt_seq
  P01053_PF00280_syns$characters <- nchar(P01053_PF00280_syns$diff)
  
  # Filter the data frame for synonymous mutants
  P01053_PF00280_syns <- P01053_PF00280_syns[P01053_PF00280_syns$characters < 12, ]
  
  
  
  # Return the results
  return(list(P01053_PF00280_syns = P01053_PF00280_syns, 
              P01053_PF00280_variants = P01053_PF00280_variants))
}

