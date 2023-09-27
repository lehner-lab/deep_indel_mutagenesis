
## 1)
## download the deep_indel_mutagenesis/Functions folder containing all the scripts with functions needed to run the analysis

# set folder location; replace ".." with your folder path
folder_location <- "~/../Functions"
script_files <- list.files(path = folder_location, pattern = "\\.R$", full.names = TRUE)

for (script_file in script_files) {
  source(script_file)
}

## 2)
## load data necessary to run all the functions in the next scripts: additional_dfs.rds

# load the additional_dfs.rds; replace ".." with you folder path
tmp=readRDS(file="~/../additional_dfs.rds")

CADD_subs = tmp[[1]]
CADD_insA = tmp[[2]]
CADD_dels = tmp[[3]]
coefficients_deletions = tmp[[4]]
coefficients_insertions = tmp[[5]]
ddmut_prediction_mean = tmp[[6]]
ddmut_predictions = tmp[[7]]
esm1v_predictions = tmp[[8]]
file_list_deletions_models = tmp[[9]]
file_list_insertions_models = tmp[[10]]
provean_dels = tmp[[11]]
provean_insA = tmp[[12]]
provean_subs = tmp[[13]]
structure_info = tmp[[14]]
all_doms_rsasa = tmp[[15]]
tsuboyama_rsasa = tmp[[16]]
validated_domains_in_vitro_ddGs_and_growthrates = tmp[[17]]
cor_tsuboyamaVStopolska = tmp[[18]]

## 3)
## download the DiMSum raw data folder set the path; replace ".." with your folder path
dimsum_file_path <- "~/../DiMSum/"

## if you run PART 1-3 (001_split_data): set folder paths to where you want to save the final processed dfs
## if you run PART4 (001_split_data): set folder paths to where you have downloaded the .rds files
# replace ".." with your folder path

output_location_aPCA <-"~/../"
output_location_bPCA <-"~/../"
output_location_colorscale <-"~/../"
output_location_tsuboyama_processed <-"~/../"

## folder path where your downloaded Tsuboyama raw data + folder with pdb_files
# replace ".." with your folder path
data_folder_tsuboyama <- "~/../Rocklin/Processed_K50_dG_datasets/"
data_folder_tsyboyama_PDBs <- "~/../AlphaFold_model_PDBs"

## folder path to save the df for running the predictor/folder path where the downloaded .rds is
output_location_ddGencoded<-"~/../"




