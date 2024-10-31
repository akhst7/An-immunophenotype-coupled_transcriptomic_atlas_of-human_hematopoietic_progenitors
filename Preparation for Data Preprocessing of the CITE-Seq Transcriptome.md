# Preparation for CITE-Seq and Single Cell Transcriptome Data Preprocessing
## Importing H5 Data and Creating a Seurat Objects

The first thing first is to download a data set from a NCBI GEO data archive under ```GSE245108```.  All the count matrix including CITE-Seq and scRNAseq matrix are stored as a H5 format file.  

![image](https://github.com/user-attachments/assets/3cc504da-e8a5-434a-8308-cc6be1f310f3)

There are dozens of h5 files consisting  titration, validation of CITE-seq antibodies and scRNAseq only data.  While in the publication, authores run analyses off a entirely merged data, for the practical purposes, the following analysis will use only several of those data from 4 study participants,  including scRNAseq of CD34hi (human HSC stem cells) and CD34+CD271+ (human BM stroma cells) sorted cells from 4 in WM_34, BM_27, WF_26, and BF_21.  
I downloaded relevant H5 files, and import them into R and generate the count matrix object by using a package called ```BPCells``` (https://bnprks.github.io/BPCells/).   
In R, there are a couple of ways to check a list of downloaded files, ```list.files``` and ```dir_ls```.  ```list.files``` ouputs names of files in the current directory whereas ```dir_ls``` shows the path of files as follows;
```
> fs::dir_ls("../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/", regexp = pattern) 
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/41590_2024_1782_MOESM2_ESM.xlsx
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_BF21-CD271_filtered_feature_bc_matrix.h5
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_BF21-CD34_filtered_feature_bc_matrix.h5
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_BF21-TNC_filtered_feature_bc_matrix.h5
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_BM27-CD271_filtered_feature_bc_matrix.h5
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_BM27-CD34_filtered_feature_bc_matrix.h5
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_BM27-TNC_filtered_feature_bc_matrix.h5
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_CD34_3_antibody_features.csv.gz
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_Level3M-titrated-cell-annotation.txt
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_WF26-CD271_filtered_feature_bc_matrix.h5
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_WF26-CD34_filtered_feature_bc_matrix.h5
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_WF26-TNC_filtered_feature_bc_matrix.h5
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_WM34-CD271_filtered_feature_bc_matrix.h5
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_WM34-CD34_filtered_feature_bc_matrix.h5
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_WM34-TNC_filtered_feature_bc_matrix.h5
```
```dir_ls``` is a command from a package ```fs```, and in this case, ```dir_ls``` is invoked without loading the ```fs``` pacakge.  Now, h5 files with **CD34** fileterd matrix will be choosen for the rest of the analysis as follows;
```
pattern <- "GSE245108_\\w.\\d.-CD34.+" # perl=TRUE
fs::dir_ls("../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/", regexp = pattern) -> path
> path
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_BF21-CD34_filtered_feature_bc_matrix.h5
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_BM27-CD34_filtered_feature_bc_matrix.h5
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_WF26-CD34_filtered_feature_bc_matrix.h5
../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/GSE245108_WM34-CD34_filtered_feature_bc_matrix.h5
```
