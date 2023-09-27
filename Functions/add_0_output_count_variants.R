### this part of the script extracts the variants with 0 output counts >5 output counts which we concider to be deleterious indel variants.
add_0_output_count_variants <- function(input_file_path) {

# Load the RData file
load(input_file_path)
  
##extract variants w. zero output counts form the Variant Data Merge File
missing_variants<-variant_data_merge %>% 
  filter(output1_e1_s1_b1_count==0 & 
           output2_e2_s1_b1_count==0 & 
           output3_e3_s1_b1_count==0)

##extract variants w. >5  input counts form the missing_variants df
missing_variants<-missing_variants %>% 
  filter(input1_e1_s0_bNA_count>5 & 
           input2_e2_s0_bNA_count>5 & 
           input3_e3_s0_bNA_count>5)

## we only find missing indel variants for: P61024_PF01111

##P61024_PF01111: find the positions at which the indels were made using wt as reference
missing_variants$P61024_PF01111=mapply(function(x, y) 
  grepRaw(diag(attr(adist(x, y, counts = TRUE), "trafos")), pattern = "S|D|I", all = TRUE),
  "IYYSDKYDDEEFEYRHVMLPKDIAKLVPKTHLMSESEWRNLGVQQSQGWVHYMIHEPEPHILLFRRP",
  missing_variants$aa_seq)

missing_variants$characters<-nchar(missing_variants$P61024_PF01111)

P61024_PF01111_missing<-missing_variants[missing_variants$characters <6,]
P61024_PF01111_missing$scaled_fitness<- 0

return(P61024_PF01111_missing)
}


