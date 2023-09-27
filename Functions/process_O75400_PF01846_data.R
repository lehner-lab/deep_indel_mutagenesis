
### this script defines the function to call the different mutational variants for O75400_PF01846
process_O75400_PF01846_data <- function(input_file_path) {
  # Load the RData file
  load(input_file_path)
  
  # Create a copy of the data for the domain
  O75400_PF01846 <- all_variants
  
  # Identify positions where sequences differ from wt_seq at the aa level
  O75400_PF01846$diff <- mapply(
    function(x, y) grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "EAKQAFKELLKEKRVPSNASWEQAMKMIINDPRYSALAKLSEKKQAFNAY",
    O75400_PF01846$aa_seq
  )
  
  # Add length of the sequence to the data frame
  O75400_PF01846$lenght <- str_length(O75400_PF01846$aa_seq)
  
  # Count the number of characters where sequences differ from wt_seq
  O75400_PF01846$characters <- nchar(O75400_PF01846$diff)
  
  # Isolate the correct sequences
  O75400_PF01846_variants <- O75400_PF01846[O75400_PF01846$characters < 14, ]
  
  # Find correct synonymous mutants for this wt
  O75400_PF01846_syns <- synonymous
  
  # Identify positions where sequences differ from wt_seq at the nt level
  O75400_PF01846_syns$diff <- mapply(
    function(x, y) grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
    "gaagcaaaacaagcatttaaagaattgttgaaagaaaaaagagttccatctaatgcatcttgggaacaagcaatgaaaatgattattaatgatccaagatattctgcattggcaaaattgtctgaaaaaaaacaagcatttaatgcatat",
    O75400_PF01846_syns$nt_seq
  )
  
  # Count the number of characters where sequences differ from wt_seq
  O75400_PF01846_syns$characters <- nchar(O75400_PF01846_syns$diff)
  
  # Filter the data frame for synonymous mutants
  O75400_PF01846_syns <- O75400_PF01846_syns[O75400_PF01846_syns$characters < 12, ]
  
 
 
  # Return the results
  return(list(O75400_PF01846_syns = O75400_PF01846_syns, 
              O75400_PF01846_variants = O75400_PF01846_variants))
}


