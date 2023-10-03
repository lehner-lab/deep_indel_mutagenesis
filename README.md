# deep_indel_mutagenesis

Welcome to the GitHub repository for: Deep indel mutagenesis reveals the impact of insertions and deletions on protein stability and function. 

# Table of contents

* **1. [Required Software](#required-software)**
* **2. [Required Data](#required-data)**
* **3. [Installation Instructions](#installation-instructions)**
* **4. [Usage](#usage)**

# Required Software

To run the deep_indel_mutagenesis pipeline you will need the following software and associated packages:

* **[_R_](https://www.r-project.org/)** (dplyr, modeest, stringr, strex, data.table, bio3d, seqinr, ggplot2, ggridges, GGally, gridExtra, viridis, writexl, matrixStats, vecsets, grpreg, glmnet, glinternet)
* **[Functions](https://github.com/lehner-lab/deep_indel_mutagenesis/tree/main/Functions)** to run the pipeline you will need to download the "Functions" folder, containing custom functions made to process the data for deep_indel_mutagenesis pipeline


# Required Data

* i) Predictions from CADD, DDMut, ESM1v, PROVEAN, ii) STRIDE and rSASA information, together with iii)  additional dfs neccesery to run the pipleline should be downloaded as "additional_dfs.rds" from  **[here](https://zenodo.org/record/8402597)** 

* DiMSum files neccesery for running the deep_indel_mutagenesis pipeline should be downloaded from **[here](https://zenodo.org/record/8402597)** 

* **[_Tsuboyama et al. 2023_](https://www.nature.com/articles/s41586-023-06328-6)**  raw data ("Tsuboyama2023_Dataset2_Dataset3_20230416.csv") and the pdb files ("AlphaFold_model_PDBs") should be downloaded **[here](https://zenodo.org/record/7992926)** 

* Pre-processed data for reproducing the figures can also be downloaded from **[here](https://zenodo.org/record/8402597)**

# Installation Instructions

Make sure you have git and conda installed and then run (expected install time <5min):

```
# Install dependencies (preferably in a fresh conda environment)
conda install -c conda-forge r-base r-dplyr r-modeest r-stringr r-strex r-data.table r-bio3d r-seqinr r-ggplot2 r-ggridges r-GGally r-gridExtra r-viridis r-writexl r-matrixStats r-vecsets r-grpreg r-glmnet r-glinternet 
```

# Usage

* **[000_load_functions](https://github.com/lehner-lab/deep_indel_mutagenesis/blob/main/000_load_functions.R)** 
  In stage 00 of the pipeline, we load and set folder location for the Required Data and load the Required Software

* **[01_split_data](https://github.com/lehner-lab/deep_indel_mutagenesis/blob/main/001_split_data.R)**
  In stage 01 of the pipeline, we process the raw DiMSum files and call the indel and substitution variants.    Furthermore we process the **[_Tsuboyama et al. 2023_](https://www.nature.com/articles/s41586-023-06328-6)**  data set for further analysis. In this script you have an option to either process the data yourself (PART1, PART2 and PART3) or directly load the processed data frames for further analysis (skip to PART4 and download **[pre-processed data](https://zenodo.org/record/8402597)**).

* **[002_figure1_main](https://github.com/lehner-lab/deep_indel_mutagenesis/blob/main/002_figure1_main.R)** Reproduce Fig. 1
* **[002_figure1_extended](https://github.com/lehner-lab/deep_indel_mutagenesis/blob/main/002_figure1_extended.R)** Reproduce Extended Fig. 1
* **[003_figure2_main](https://github.com/lehner-lab/deep_indel_mutagenesis/blob/main/003_figure2_main.R)** Reproduce Fig. 2
* **[003_figure2_extended](https://github.com/lehner-lab/deep_indel_mutagenesis/blob/main/003_figure2_extended.R)** Reproduce Extended Fig. 2
* **[004_figure3_main](https://github.com/lehner-lab/deep_indel_mutagenesis/blob/main/004_figure3_main.R)** Reproduce Fig. 3
* **[004_figure3_extended](https://github.com/lehner-lab/deep_indel_mutagenesis/blob/main/004_figure3_main.R)** Reproduce Extended Fig. 3
* **[005_figure4_main](https://github.com/lehner-lab/deep_indel_mutagenesis/blob/main/005_figure4_main.R)** Reproduce Fig. 4
* **[005_figure4_extended](https://github.com/lehner-lab/deep_indel_mutagenesis/blob/main/005_figure4_main.R)** Reproduce Extended Fig. 4
* **[006_figure5_main](https://github.com/lehner-lab/deep_indel_mutagenesis/blob/main/006_figure5_main.R)** Reproduce Fig. 5
* **[006_figure5_extended](https://github.com/lehner-lab/deep_indel_mutagenesis/blob/main/006_figure5_extended.R)** Reproduce Extended Fig. 5

000_load_functions and 01_split_data should be run first if i) you want to reproduce the final data frames from raw data and ii) if you want to work with the pre-processed data frames. 


Scripts to reproduce the figures can be run independently as long as they are run in pairs (if you want to run 003_figure2_extended, you should run 003_figure2_main first). 





  



