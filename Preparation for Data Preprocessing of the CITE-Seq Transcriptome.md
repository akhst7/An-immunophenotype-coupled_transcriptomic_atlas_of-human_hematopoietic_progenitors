# Preparation for CITE-Seq and Single Cell Transcriptome Data Preprocessing
## Importing H5 Data and Creating a Seurat Objects

The first thing first is to download a data set from a NCBI GEO data archive under ```GSE245108```.  All the count matrix including CITE-Seq and scRNAseq matrix are stored as a H5 format file.  

![image](https://github.com/user-attachments/assets/3cc504da-e8a5-434a-8308-cc6be1f310f3)

There are dozens of h5 files consisting  titration, validation of CITE-seq antibodies and scRNAseq only data.  While in the publication, authores run analyses off a entirely merged data, for the practical purposes, the following analysis will use only several of those data from 4 study participants,  including scRNAseq of CD34hi (human HSC stem cells) and CD34+CD271+ (human BM stroma cells) sorted cells from 4 in WM_34, BM_27, WF_26, and BF_21.  
I downloaded relevant H5 files, and import them into R and generate the count matrix object by using a package called ```BPCells``` (https://bnprks.github.io/BPCells/).   
