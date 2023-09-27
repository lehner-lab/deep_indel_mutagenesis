
### this script defines the function to call the different mutational variants for P32081_PF00313
process_P32081_PF00313_data <- function(input_file_path) {
  # Load the RData file
  load(input_file_path)
  
  # Create a copy of the data for the domain
  P32081_PF00313 <- all_variants
  
  # Identify positions where sequences differ from wt_seq at the aa level
  P32081_PF00313$diff <- mapply(
    function(x, y) grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "EGKVKWFNSEKGFGFIEVEGQDDVFVHFSAIQGEGFKTLEEGQAVSFEIVEGNRGPQAANVTK",
    P32081_PF00313$aa_seq
  )
  
  # Add length of the sequence to the data frame
  P32081_PF00313$lenght <- str_length(P32081_PF00313$aa_seq)
  
  # Count the number of characters where sequences differ from wt_seq
  P32081_PF00313$characters <- nchar(P32081_PF00313$diff)
  
  # Isolate the correct sequences
  P32081_PF00313_variants <- P32081_PF00313[P32081_PF00313$characters < 14, ]
  
  # Find correct synonymous mutants for this wt
  P32081_PF00313_syns <- synonymous
  
  # Identify positions where sequences differ from wt_seq at the nt level
  P32081_PF00313_syns$diff <- mapply(
    function(x, y) grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "gaaggtaaagttaaatggtttaattctgaaaaaggttttggttttattgaagttgaaggtcaagatgatgtttttgttcatttttctgcaattcaaggtgaaggttttaaaactttggaagaaggtcaagcagtttcttttgaaattgttgaaggtaatagaggtccacaagcagcaaatgttactaaa",
    P32081_PF00313_syns$nt_seq
  )
  
  # Count the number of characters where sequences differ from wt_seq
  P32081_PF00313_syns$characters <- nchar(P32081_PF00313_syns$diff)
  
  # Filter the data frame for synonymous mutants
  P32081_PF00313_syns <- P32081_PF00313_syns[P32081_PF00313_syns$characters < 12, ]
  
  
  
  # Return the results
  return(list(P32081_PF00313_syns = P32081_PF00313_syns, 
              P32081_PF00313_variants = P32081_PF00313_variants))
}

