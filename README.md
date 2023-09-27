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

* i) Predictions from CADD, DDMut, ESM1v, PROVEAN, ii) STRIDE and rSASA information, together with iii)  additional dfs neccesery to run the pipleline should be downloaded as "additional_dfs.rds" from  **[here](https://crgcnag-my.sharepoint.com/personal/mtopolska_crg_es/_layouts/15/onedrive.aspx?login_hint=mtopolska%40crg%2Ees&id=%2Fpersonal%2Fmtopolska%5Fcrg%5Fes%2FDocuments%2FTopolska%5Fetal%5Fdeep%5Findel%5Fmutagenesis&view=0)** 

* DiMSum files neccesery for running the deep_indel_mutagenesis pipeline should be downloaded from **[here](https://crgcnag-my.sharepoint.com/personal/mtopolska_crg_es/_layouts/15/onedrive.aspx?login_hint=mtopolska%40crg%2Ees&id=%2Fpersonal%2Fmtopolska%5Fcrg%5Fes%2FDocuments%2FTopolska%5Fetal%5Fdeep%5Findel%5Fmutagenesis%2FDiMSum&view=0)** 

* **[_Tsuboyama et al. 2023_](https://www.nature.com/articles/s41586-023-06328-6)**  raw data ("Tsuboyama2023_Dataset2_Dataset3_20230416.csv") and the pdb files ("AlphaFold_model_PDBs") should be downloaded **[here](https://zenodo.org/record/7992926)** 

* Pre-processed data for reproducing the figures can also be downloaded from **[here](https://zenodo.org/record/7992926](https://crgcnag-my.sharepoint.com/personal/mtopolska_crg_es/_layouts/15/onedrive.aspx?login_hint=mtopolska%40crg%2Ees&id=%2Fpersonal%2Fmtopolska%5Fcrg%5Fes%2FDocuments%2FTopolska%5Fetal%5Fdeep%5Findel%5Fmutagenesis%2Fpreprocessed%5Fdata&view=0)**

# Installation Instructions

Make sure you have git and conda installed and then run (expected install time <5min):

```
# Install dependencies (preferably in a fresh conda environment)
conda install -c conda-forge r-base r-dplyr r-modeest r-stringr r-strex r-data.table r-bio3d r-seqinr r-ggplot2 r-ggridges r-GGally r-gridExtra r-viridis r-writexl r-matrixStats r-vecsets r-grpreg r-glmnet r-glinternet 
```

# Usage

* **[00_load_functions]()** 
  In stage 00 of the pipeline, we load and set folder location for the Required Data and load the Required Software

* **[01_split_data]()**
  In stage 01 of the pipeline, we process the raw DiMSum files and call the indel and substitution variants.    Furthermore we process the **[_Tsuboyama et al. 2023_](https://www.nature.com/articles/s41586-023-06328-6)**  data set for further analysis. In this script you have an option to either process the data yourself (PART1, PART2 and PART3) or directly load the processed data frames for further analysis (skip to PART4 and download **[pre-processed data](https://crgcnag-my.sharepoint.com/personal/mtopolska_crg_es/_layouts/15/onedrive.aspx?login_hint=mtopolska%40crg%2Ees&id=%2Fpersonal%2Fmtopolska%5Fcrg%5Fes%2FDocuments%2FTopolska%5Fetal%5Fdeep%5Findel%5Fmutagenesis%2Fpreprocessed%5Fdata&view=0)**).

* **[002_figure1_main]()** Reproduce Fig. 1
* **[002_figure1_extended]()** Reproduce Extended Fig. 1
* **[003_figure2_main]()** Reproduce Fig. 2
* **[003_figure2_extended]()** Reproduce Extended Fig. 2
* **[004_figure3_main]()** Reproduce Fig. 3
* **[004_figure3_extended]()** Reproduce Extended Fig. 3
* **[005_figure4_main]()** Reproduce Fig. 4
* **[005_figure4_extended]()** Reproduce Extended Fig. 4
* **[006_figure5_main]()** Reproduce Fig. 5
* **[006_figure5_extended]()** Reproduce Extended Fig. 5

Stage 00 and 01 should be run first. The following stages can be ran independently as long as they are ran in pairs (if you want to run 003_figure2_extended, you should run 003_figure2_main first). 





  



