
## 1)
## download the all_functions folder containing all the scripts with functions needed to run the analysis

folder_location <- "~/Desktop/mount1/indel2fold/many_domains/joined_analysis/doms_joined_analysis_fold_longer_selection/analysis/code_manuscript/functions/all_functions"
script_files <- list.files(path = folder_location, pattern = "\\.R$", full.names = TRUE)

for (script_file in script_files) {
  source(script_file)
}

## 2)
## load data necessary to run all the functions in the next scripts: additional_dfs.rds
tmp=readRDS(file="~/Desktop/mount1/indel2fold/many_domains/joined_analysis/doms_joined_analysis_fold_longer_selection/analysis/code_manuscript/functions/additional_dfs.rds")

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
### download the analysis folder and set the different paths. 
## folder where the dimsum files are
dimsum_file_path <- "~/Desktop/mount1/indel2fold/many_domains/joined_analysis/doms_joined_analysis_fold_longer_selection/analysis/code_manuscript/functions/DimSum/"

## folder where you want to save the final data frames
output_location_aPCA<-"~/Desktop/mount1/indel2fold/many_domains/joined_analysis/doms_joined_analysis_fold_longer_selection/analysis/code_manuscript/functions/"
output_location_bPCA<-"~/Desktop/mount1/indel2fold/many_domains/joined_analysis/doms_joined_analysis_fold_longer_selection/analysis/code_manuscript/functions/"
output_location_colorscale<-"~/Desktop/mount1/indel2fold/many_domains/joined_analysis/doms_joined_analysis_fold_longer_selection/analysis/code_manuscript/functions/"
output_location_tsuboyama_processed<-"~/Desktop/mount1/indel2fold/many_domains/joined_analysis/doms_joined_analysis_fold_longer_selection/analysis/code_manuscript/functions/"

## folder where your downloaded Tsuboyama data + pdb files are
data_folder_tsuboyama<-"~/Desktop/indel_pub/External_data/Rocklin/Processed_K50_dG_datasets/"
data_folder_tsyboyama_PDBs<-"~/Desktop/mount1/Rocklin/Rocklin_pdbs/"

## where to save the data encoded to run the predictor
output_location_ddGencoded<-"~/Desktop/mount1/indel2fold/many_domains/joined_analysis/doms_joined_analysis_fold_longer_selection/analysis/code_manuscript/functions/"




