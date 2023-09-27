## load neccessery packages

library (dplyr)
library(modeest)
library(stringr)
library(strex)
library(data.table)
library(bio3d)
library(seqinr)


# PART 1: This script processes the raw aPCA data from the 9 domain scan. The final output is the "scaled_variants_aPCA" df.
# PART 2: This script processes the raw bPCA data from the 2 domain scan. The final output is the "scaled_variants_bPCA" df.
# PART 3: This script processes the Tsuboyama data from: Tsuboyama, K., Dauparas, J., Chen, J. et al. Mega-scale experimental analysis of protein folding stability in biology and design. Nature 620, 434–444 (2023). https://doi.org/10.1038/s41586-023-06328-6

# PART 4: skip to PART 4 if you want to work with the pre-processed data from PART1+2+3

##########################################
#############   PART 1  ##################
##########################################

## this part of the script is used to call variants and generate a df for data analysis of the aPCA 9 domain dataset.

# Call the function to import the "missing variants" from the Variant Data Merge file; variants with 0 counts in output but >6 counts in the input
# we concider those as totally deleterious indels. 
input_file_path <- paste(dimsum_file_path, "aPCA_domains_variant_data_merge.RData",sep="/")
P61024_PF01111_missing <- add_0_output_count_variants(input_file_path)


## 1)
# Call the function to process data for O75400_PF01846 (FBP11-FF1) using the Dimsum raw data from the 7 domain selection
# here all mutation variants are called
input_file_path <- paste(dimsum_file_path, "aPCA_domains_fitness_replicates.RData",sep="/")
results <- process_O75400_PF01846_data(input_file_path)
O75400_PF01846_syns=results[[1]]
O75400_PF01846_variants=results[[2]]


# Call the function to normalise the data 
results <- normalise_O75400_PF01846_data()

O75400_PF01846_allvariants = results[[1]]
O75400_PF01846_syns = results[[2]]
O75400_PF01846_CX = results[[3]]
O75400_PF01846_CXX = results[[4]]
O75400_PF01846_CXXX = results[[5]]
O75400_PF01846_singleDEL = results[[6]]
O75400_PF01846_doubleDEL = results[[7]]
O75400_PF01846_tripleDEL = results[[8]]
O75400_PF01846_Delsub = results[[9]]
O75400_PF01846_Delsub_2 = results[[10]]


## 2)
# Call the function to process data for P0A9X9_PF00313 (CSPA-CSD) using the Dimsum raw data from the 7 domain selection
# here all mutation variants are called
input_file_path <- paste(dimsum_file_path, "aPCA_domains_fitness_replicates.RData",sep="/")
results <- process_P0A9X9_PF00313_data(input_file_path)
P0A9X9_PF00313_syns=results[[1]]
P0A9X9_PF00313_variants=results[[2]]


# Call the function to normalise the data 
results <- normalise_P0A9X9_PF00313_data()

P0A9X9_PF00313_allvariants = results[[1]]
P0A9X9_PF00313_syns = results[[2]]
P0A9X9_PF00313_CX = results[[3]]
P0A9X9_PF00313_CXX = results[[4]]
P0A9X9_PF00313_CXXX = results[[5]]
P0A9X9_PF00313_singleDEL = results[[6]]
P0A9X9_PF00313_doubleDEL = results[[7]]
P0A9X9_PF00313_tripleDEL = results[[8]]
P0A9X9_PF00313_Delsub = results[[9]]
P0A9X9_PF00313_Delsub_2 = results[[10]]

## 3)
# Call the function to process data for P01053_PF00280 (CI2A-PIN1) using the Dimsum raw data from the 7 domain selection
# here all mutation variants are called
input_file_path <- paste(dimsum_file_path, "aPCA_domains_fitness_replicates.RData",sep="/")
results <- process_P01053_PF00280_data(input_file_path)
P01053_PF00280_syns=results[[1]]
P01053_PF00280_variants=results[[2]]


# Call the function to normalise the data
results <- normalise_P01053_PF00280_data()

P01053_PF00280_allvariants = results[[1]]
P01053_PF00280_syns = results[[2]]
P01053_PF00280_CX = results[[3]]
P01053_PF00280_CXX = results[[4]]
P01053_PF00280_CXXX = results[[5]]
P01053_PF00280_singleDEL = results[[6]]
P01053_PF00280_doubleDEL = results[[7]]
P01053_PF00280_tripleDEL = results[[8]]
P01053_PF00280_Delsub = results[[9]]

## 4)
# Call the function to process data for P02417_PF01281 (BL17-NTL9) using the Dimsum raw data from the 7 domain selection
# here all mutation variants are called
input_file_path <- paste(dimsum_file_path, "aPCA_domains_fitness_replicates.RData",sep="/")
results <- process_P02417_PF01281_data(input_file_path)
P02417_PF01281_syns=results[[1]]
P02417_PF01281_variants=results[[2]]


# Call the function to normalise the data 
results <- normalise_P02417_PF01281_data()

P02417_PF01281_allvariants = results[[1]]
P02417_PF01281_syns = results[[2]]
P02417_PF01281_CX = results[[3]]
P02417_PF01281_CXX = results[[4]]
P02417_PF01281_CXXX = results[[5]]
P02417_PF01281_singleDEL = results[[6]]
P02417_PF01281_doubleDEL = results[[7]]
P02417_PF01281_tripleDEL = results[[8]]
P02417_PF01281_Delsub = results[[9]]
P02417_PF01281_Delsub_2 = results[[10]]

## 5)
# Call the function to process data for P02640_PF02209 (VIL1-HP) using the Dimsum raw data from the 7 domain selection
# here all mutation variants are called
input_file_path <- paste(dimsum_file_path, "aPCA_domains_fitness_replicates.RData",sep="/")
results <- process_P02640_PF02209_data(input_file_path)
P02640_PF02209_syns=results[[1]]
P02640_PF02209_variants=results[[2]]


# Call the function to normalise the data 
results <- normalise_P02640_PF02209_data()

P02640_PF02209_allvariants = results[[1]]
P02640_PF02209_syns = results[[2]]
P02640_PF02209_CX = results[[3]]
P02640_PF02209_CXX = results[[4]]
P02640_PF02209_CXXX = results[[5]]
P02640_PF02209_singleDEL = results[[6]]
P02640_PF02209_doubleDEL = results[[7]]
P02640_PF02209_tripleDEL = results[[8]]
P02640_PF02209_Delsub = results[[9]]
P02640_PF02209_Delsub_2 = results[[10]]

## 6)
# Call the function to process data for P32081_PF00313 (CSPB-CSD) using the Dimsum raw data from the 7 domain selection
# here all mutation variants are called
input_file_path <- paste(dimsum_file_path, "aPCA_domains_fitness_replicates.RData",sep="/")
results <- process_P32081_PF00313_data(input_file_path)
P32081_PF00313_syns=results[[1]]
P32081_PF00313_variants=results[[2]]


# Call the function to normalise the data 
results <- normalise_P32081_PF00313_data()

P32081_PF00313_allvariants = results[[1]]
P32081_PF00313_syns = results[[2]]
P32081_PF00313_CX = results[[3]]
P32081_PF00313_CXX = results[[4]]
P32081_PF00313_CXXX = results[[5]]
P32081_PF00313_singleDEL = results[[6]]
P32081_PF00313_doubleDEL = results[[7]]
P32081_PF00313_tripleDEL = results[[8]]
P32081_PF00313_Delsub = results[[9]]
P32081_PF00313_Delsub_2 = results[[10]]

## 7)
# Call the function to process data for P61024_PF01111 (CKS1) using the Dimsum raw data from the 7 domain selection
# here all mutation variants are called
input_file_path <- paste(dimsum_file_path, "aPCA_domains_fitness_replicates.RData",sep="/")
results <- process_P61024_PF01111_data(input_file_path)
P61024_PF01111_syns=results[[1]]
P61024_PF01111_variants=results[[2]]


# Call the function to normalise the data 
results <- normalise_P61024_PF01111_data()

P61024_PF01111_allvariants = results[[1]]
P61024_PF01111_syns = results[[2]]
P61024_PF01111_CX = results[[3]]
P61024_PF01111_CXX = results[[4]]
P61024_PF01111_CXXX = results[[5]]
P61024_PF01111_singleDEL = results[[6]]
P61024_PF01111_doubleDEL = results[[7]]
P61024_PF01111_tripleDEL = results[[8]]
P61024_PF01111_Delsub = results[[9]]
P61024_PF01111_Delsub_2 = results[[10]]


## 8)
# Call the function to process data for GRB2-SH3 using the Dimsum raw data from the aPCA selection for GRB2-SH3
# here all mutation variants are called
results <- process_grb2_fold_data()

grb2_fold_CX = results[[1]]
grb2_fold_CXX = results[[2]]
grb2_fold_CNN = results[[3]]
grb2_fold_insAA = results[[4]]
grb2_fold_insCC = results[[5]]
grb2_fold_CXXX = results[[6]]
grb2_fold_CNNN = results[[7]]
grb2_fold_singles = results[[8]]
grb2_fold_allins = results[[9]]
grb2_fold_Delsub = results[[10]]
grb2_fold_wildtype = results[[11]]
grb2_fold_synonymous = results[[12]]
grb2_fold_singleDEL = results[[13]]
grb2_fold_doubleDEL = results[[14]]
grb2_fold_tripleDEL = results[[15]]


# Call the function to normalise the data 
results <- normalise_grb2_fold_data()

grb2_fold_allvariants = results[[1]] 
grb2_fold_CX = results[[2]] 
grb2_fold_CXX = results[[3]] 
grb2_fold_CXXX = results[[4]] 
grb2_fold_CNN = results[[5]] 
grb2_fold_CNNN = results[[6]] 
grb2_fold_insCC = results[[7]] 
grb2_fold_insAA = results[[8]] 
grb2_fold_allins = results[[9]] 
grb2_fold_synonymous = results[[10]] 
grb2_fold_singles = results[[11]] 
grb2_fold_Delsub = results[[12]] 
grb2_fold_Delsub_2 = results[[13]] 
grb2_fold_singleDEL = results[[14]]
grb2_fold_doubleDEL = results[[15]]
grb2_fold_tripleDEL = results[[16]]


## 9)
# Call the function to process data for PSD95-PDZ3 using the Dimsum raw data from the aPCA selection for PSD95-PDZ3
# here all mutation variants are called
results <- process_pdz3_fold_data()

pdz3_fold_CX = results[[1]]
pdz3_fold_CXX = results[[2]]
pdz3_fold_CNN = results[[3]]
pdz3_fold_insAA = results[[4]]
pdz3_fold_insCC = results[[5]]
pdz3_fold_CXXX = results[[6]]
pdz3_fold_CNNN = results[[7]]
pdz3_fold_singles = results[[8]]
pdz3_fold_allins = results[[9]]
pdz3_fold_Delsub = results[[10]]
pdz3_fold_wildtype = results[[11]]
pdz3_fold_singleDEL = results[[12]]
pdz3_fold_doubleDEL = results[[13]]
pdz3_fold_tripleDEL = results[[14]]
pdz3_fold_synonymous = results[[15]]

# Call the function to normalise the data 
results <- normalise_pdz3_fold_data()

pdz3_fold_allvariants = results[[1]] 
pdz3_fold_CX = results[[2]] 
pdz3_fold_CXX = results[[3]] 
pdz3_fold_CXXX = results[[4]] 
pdz3_fold_CNN = results[[5]] 
pdz3_fold_CNNN = results[[6]] 
pdz3_fold_insCC = results[[7]] 
pdz3_fold_insAA = results[[8]] 
pdz3_fold_allins = results[[9]] 
pdz3_fold_synonymous = results[[10]] 
pdz3_fold_singles = results[[11]] 
pdz3_fold_Delsub = results[[12]] 
pdz3_fold_Delsub_2 = results[[13]] 
pdz3_fold_singleDEL = results[[14]] 
pdz3_fold_doubleDEL = results[[15]] 
pdz3_fold_tripleDEL = results[[16]] 

## This function adds domain name to all aPCA dfs and joins all dfs into a common df. 
results <- add_domain_name_and_join_dfs()
scaled_variants_aPCA <- results[[1]] 
scaled_variants_aPCA$Pos<-as.numeric(scaled_variants_aPCA$Pos)

###  structural info obtained by STRIDE: http://webclu.bio.wzw.tum.de/cgi-bin/stride/stridecgi.py
# loaded in 00_load_function script as: structure_info

##  rsasa information obtained with Pymol
# loaded in 00_load_function script as: all_doms_rsasa


## merge to add structural (stride) and rSASA (pymol) info
scaled_variants_aPCA <- merge(scaled_variants_aPCA,
                            structure_info[,c("WT_aa", "Pos", "Structure_name", "domain")],
                            by=c("Pos", "domain"),
                            all.x = T)

scaled_variants_aPCA <- merge(scaled_variants_aPCA,
                              all_doms_rsasa[,c("Pos", "rSASA", "domain")],
                              by=c("Pos", "domain"),
                              all.x = T)

## define location where you want to save this .rds
output_location<-output_location_aPCA
output_name<-"scaled_variants_aPCA.rds"

saveRDS(list(scaled_variants_aPCA), 
        file=paste(output_location,
                   output_name, sep="/"))

##########################################
#############   PART 2  ##################
##########################################
## this script is used to call variants and generate a df for data analysis of the bPCA for 2 domain dataset.

## 1)
# Call the function to process data for GRB2-SH3 using the Dimsum raw data from the bPCA for GRB2-SH3
# here all mutation variants are called
results <- process_grb2_bind_data()

grb2_bind_CX = results[[1]]
grb2_bind_CXX = results[[2]]
grb2_bind_CNN = results[[3]]
grb2_bind_insAA = results[[4]]
grb2_bind_insCC = results[[5]]
grb2_bind_CXXX = results[[6]]
grb2_bind_CNNN = results[[7]]
grb2_bind_singles = results[[8]]
grb2_bind_allins = results[[9]]
grb2_bind_Delsub = results[[10]]
grb2_bind_wildtype = results[[11]]
grb2_bind_synonymous = results[[12]]
grb2_bind_singleDEL = results[[13]]
grb2_bind_doubleDEL = results[[14]]
grb2_bind_tripleDEL = results[[15]]


# Call the function to normalise the data 
results <- normalise_grb2_bind_data()

grb2_bind_allvariants = results[[1]] 
grb2_bind_CX = results[[2]] 
grb2_bind_CXX = results[[3]] 
grb2_bind_CXXX = results[[4]] 
grb2_bind_CNN = results[[5]] 
grb2_bind_CNNN = results[[6]] 
grb2_bind_insCC = results[[7]] 
grb2_bind_insAA = results[[8]] 
grb2_bind_allins = results[[9]] 
grb2_bind_synonymous = results[[10]] 
grb2_bind_singles = results[[11]] 
grb2_bind_Delsub = results[[12]] 
grb2_bind_Delsub_2 = results[[13]] 
grb2_bind_singleDEL = results[[14]]
grb2_bind_doubleDEL = results[[15]]
grb2_bind_tripleDEL = results[[16]]


## 2)
# Call the function to process data for PSD95-PDZ3 using the Dimsum raw data from the bPCA for PSD95-PDZ3
# here all mutation variants are called
results <- process_pdz3_bind_data()

pdz3_bind_CX = results[[1]]
pdz3_bind_CXX = results[[2]]
pdz3_bind_CNN = results[[3]]
pdz3_bind_insAA = results[[4]]
pdz3_bind_insCC = results[[5]]
pdz3_bind_CXXX = results[[6]]
pdz3_bind_CNNN = results[[7]]
pdz3_bind_singles = results[[8]]
pdz3_bind_allins = results[[9]]
pdz3_bind_Delsub = results[[10]]
pdz3_bind_wildtype = results[[11]]
pdz3_bind_singleDEL = results[[12]]
pdz3_bind_doubleDEL = results[[13]]
pdz3_bind_tripleDEL = results[[14]]
pdz3_bind_synonymous = results[[15]]

# Call the function to normalise the data 
results <- normalise_pdz3_bind_data()

pdz3_bind_allvariants = results[[1]] 
pdz3_bind_CX = results[[2]] 
pdz3_bind_CXX = results[[3]] 
pdz3_bind_CXXX = results[[4]] 
pdz3_bind_CNN = results[[5]] 
pdz3_bind_CNNN = results[[6]] 
pdz3_bind_insCC = results[[7]] 
pdz3_bind_insAA = results[[8]] 
pdz3_bind_allins = results[[9]] 
pdz3_bind_synonymous = results[[10]] 
pdz3_bind_singles = results[[11]] 
pdz3_bind_Delsub = results[[12]] 
pdz3_bind_Delsub_2 = results[[13]] 
pdz3_bind_singleDEL = results[[14]] 
pdz3_bind_doubleDEL = results[[15]] 
pdz3_bind_tripleDEL = results[[16]] 

## This function adds domain name to all aPCA dfs and joins all dfs into a common df. 
results <- add_domain_name_and_join_dfs_bPCA()
scaled_variants_bPCA <- results[[1]] 
scaled_variants_bPCA$Pos<-as.numeric(scaled_variants_bPCA$Pos)

## merge to add structural (stride) and rSASA (pymol) info
scaled_variants_bPCA <- merge(scaled_variants_bPCA,
                              structure_info[,c("WT_aa", "Pos", "Structure_name", "domain")],
                              by=c("Pos", "domain"),
                              all.x = T)

scaled_variants_bPCA <- merge(scaled_variants_bPCA,
                              all_doms_rsasa[,c("Pos", "rSASA", "domain")],
                              by=c("Pos", "domain"),
                              all.x = T)

## define location where you want to save this .rds
output_location<-output_location_bPCA
output_name<-"scaled_variants_bPCA.rds"

saveRDS(list(scaled_variants_bPCA), 
        file=paste(output_location,
                   output_name, sep="/"))

## download the joined df to set common color scale for heatmaps. 
color_scale <- c(scaled_variants_aPCA$scaled_fitness,
                 scaled_variants_bPCA$scaled_fitness)

## define location where you want to save this .rds
output_location<-output_location_colorscale
output_name<-"color_scale.rds"

saveRDS(list(color_scale), 
        file=paste(output_location,
                   output_name, sep="/"))


##########################################
#############   PART 3  ##################
##########################################

## this script is used to download and process Tsuboyama et al., dataset for further analysis.
## to use this script you need to download a) Tsuboyama2023_Dataset2_Dataset3_20230416.csv and 
#                                          b) AlphaFold_model_PDBs from: https://zenodo.org/record/7992926 and 
#                                          c) the rSASA information obtained from the AlphaFold_model_PDBs using PyMol: tsuboyama_rsasa.rds (loaded in 00_load_functions step)

# set folder where your Tsuboyama2023_Dataset2_Dataset3_20230416.csv is located
data_folder <- paste(data_folder_tsuboyama, "Tsuboyama2023_Dataset2_Dataset3_20230416.csv",sep="")

## script to isolate natural domains into a list split based on pdbname (~ this step can take a couple of minutes to run)
result <- isolate_natural_doms(data_folder)

data_list_nat_doms = result[[1]]
pdbnames = result[[2]]

## run stride for all Tsuboyama pdb domains. (~ this step can take a couple of minutes to run)
# set folder where your AlphaFold_model_PDBs is located
data_folder <- paste(data_folder_tsyboyama_PDBs, "AlphaFold_model_PDBs",sep="")

result <- run_stride(data_folder)
strucutre_tsuboyama_list = result

## merge stride and pdb information together with data. 
# extract pdbnames
names_nat_doms <- names(data_list_nat_doms)
names_struc_info <- names(strucutre_tsuboyama_list)

list_tsuboyama_nat_doms_strucutre<-c()
for (i in names_nat_doms) {
  df1 <- data_list_nat_doms[[i]]
  df2 <- strucutre_tsuboyama_list[[i]]
  result <- merge(df1, df2, by = c("Pos"))
  list_tsuboyama_nat_doms_strucutre[[i]] <- result
}

## create final df for analysis
result <- final_tsuboyama_df()

tsuboyama_nat_doms_all <- result[[1]]
tsuboyama_nat_doms_structure_subs <- result[[2]]

## add structural feature information to both dfs. (~ this step can take a couple of minutes to run)
##  rsasa information obtained with Pymol
# loaded in 00_load_function script: tsuboyama_rsasa
result <- structural_features_tsuboyama()

tsuboyama_sec_struc <- result[[1]]

## merge final df. 
tsuboyama_nat_doms_all <- merge(tsuboyama_nat_doms_all,
                                tsuboyama_sec_struc[,-c(3:4)],
                                by = c("Pos", "pdb_name"))

tsuboyama_nat_doms_structure_subs <- merge(tsuboyama_nat_doms_structure_subs,
                                           tsuboyama_sec_struc[,-c(3:4)],
                                           by = c("Pos", "pdb_name"))

## define location where you want to save this .rds
output_location <- output_location_tsuboyama_processed
output_name<-"tsuboyama_nat_doms_all.rds"

saveRDS(list(tsuboyama_nat_doms_all,
             tsuboyama_nat_doms_structure_subs), 
        file=paste(output_location,
                   output_name, sep="/"))










####################################################################################
########### skip to here if you want to work with the pre-processed data. ##########
####################################################################################

library (dplyr)
library(modeest)
library(stringr)
library(strex)
library(data.table)
library(bio3d)
library(seqinr)

##########################################
#############   PART 4  ##################
##########################################

## download data for aPCA of 9 domains, filename: scaled_variants_aPCA.rds
# define the location where the scaled_variants_aPCA.rds is downloaded and saved
output_location <- output_location_aPCA
output_name<-"scaled_variants_aPCA.rds"
# read the df
tmp<-readRDS(file=paste(output_location,
                   output_name, sep="/"))

scaled_variants_aPCA=tmp[[1]]


## download data for bPCA of 2 domains, filename: scaled_variants_bPCA
# define the location where the scaled_variants_aPCA.rds is downloaded and saved
output_location <- output_location_bPCA
output_name<-"scaled_variants_bPCA.rds"
# read the df
tmp<-readRDS(file=paste(output_location,
                        output_name, sep="/"))

scaled_variants_bPCA=tmp[[1]]


## download data for coomon color scale across all heatmaps, filename: color_scale
# define the location where the color_scale.rds is downloaded and saved
output_location <- output_location_colorscale
output_name<-"color_scale.rds"
# read the df
tmp<-readRDS(file=paste(output_location,
                        output_name, sep="/"))

color_scale=tmp[[1]]

## download processed data for Tsuboyama et al. mega-dms, filename: scaled_variants_aPCA.rds
# define the location where the scaled_variants_aPCA.rds is downloaded and saved
output_location <- output_location_tsuboyama_processed
output_name<-"tsuboyama_nat_doms_all.rds"
# read the df
tmp<-readRDS(file=paste(output_location,
                        output_name, sep="/"))

tsuboyama_nat_doms_all = tmp[[1]]
tsuboyama_nat_doms_structure_subs = tmp[[2]]

rm(tmp)






