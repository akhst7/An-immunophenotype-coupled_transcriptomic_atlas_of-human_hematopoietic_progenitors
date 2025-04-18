# CITE-Seq and Single Cell Transcriptome Data Preprocessing
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
These are h5 matrix files for sorted CD34+ cells from human BM.  These files contain both GEX(gene expression) and Antibody(CITE-seq) matrix and for creating Seurat objs, the first the GEX matrix needs to be imported.  Then, the CITE-seq matrix can be added into the Seurat obj as an new assay.  Although individual Seurat objs could be created by repeating a following with changing path[x] from 1 to for 4 times and later merge them into the one giant Seurat obj, perhaps the simpler way to accomplish the same is to use ```lapply``` looping over a list of ```path``` as follows;
```
GEX_su<-lapply(1:4, function(x){
dir=path[x]
tmp<-Read10X_h5(dir)
CreateSeuratObject(tmp$`Gene Expression`)
  }
)
```
Names of individual items of ```GEX_su``` should be added to make the rest of lives easier  by;
```
names(GEX_su)<-c("BF21", "BM27", "WF26", "WM34")
```
Then each imported h5 matrix in the ```GEX_SU``` will be merged by;
```
merge.su<-merge(GEX_su[[1]], GEX_su[2:4])
```
In the SeuratV5 obj, these merged count matrix will be added as ```layers```, and in fact merge.su has 4 groups with over 6000 cells.  
```
Formal class 'Assay5' [package "SeuratObject"] with 8 slots
  ..@ layers    :List of 9
  .. ..$ counts.BF12_CD34:Formal class 'MatrixSubset' [package "BPCells"] with 7 slots
  .. ..$ counts.BM27_CD34:Formal class 'MatrixSubset' [package "BPCells"] with 7 slots
  .. ..$ counts.WF26_CD34:Formal class 'MatrixSubset' [package "BPCells"] with 7 slots
  .. ..$ counts.WM34_CD34:Formal class 'MatrixSubset' [package "BPCells"] with 7 slots

> table(merge.su$group)

 BF21  BM27  WF26  WM34 
 6622  8209 11698  8213 

```
Overall, this Seurat object is **34742** total number of cells, and this should makes this Seurat object relatively large. 
```
format(object.size(merge.su), units = "GB")
[1] "5.1 Gb"
```
If using R in the cloud, whether using an institutional HPC or Cloud services such as AWS  accessible with an a large amount of RAM, having large Seurat objects wont matter that much however, local computer with limited resources (e.g. CPU and RAM), will hit the momory limit rather quickly. The magic to get around this is to use a package, called ```BPCells``` (https://github.com/bnprks/BPCells).  A following is a short description of ```BPCells``` functions;
```
BPCells is a package for high performance single cell analysis on RNA-seq and ATAC-seq datasets. It can analyze a 1.3M cell dataset with 2GB of RAM in under 10 minutes. This makes analysis of million-cell datasets practical on a laptop.

BPCells provides:

Efficient storage of single cell datasets via bitpacking compression
Fast, disk-backed RNA-seq and ATAC-seq data processing powered by C++
Downstream analysis such as marker genes, and clustering
Interoperability with AnnData, 10x datasets, R sparse matrices, and GRanges
Additionally, BPCells exposes its optimized data processing infrastructure for use in scaling 3rd party single cell tools (e.g. Seurat)
```
A Seurat web site has a tutorial showing how to deploy a BPCell matrix in their workflow (https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette).  Basically, BPCells allows to leave the large portion of data in the disk with **bitpacking compression**  rather than importing it to the RAM, thereby saving a memory space.  Bitpacking compression is written C++, which makes the computational process much faster.   ```BPCells``` is not compatible with a Bioconductor SingleCellExperiment platform, however there is a sister packages, ```BPCellsArray``` which works for the SingleCellExperiment objects.  
A following section describes how to utilize ```BPCells``` for the Seurat objects.

## Implementing BPCells in Seurat 
For a multiple samples, the easiest way to incorporate BPCells matrix is to generate individual Surate objs with BPCell count matrix and merge them later.  In this example, there are 4 samples (there are total of 13 samples) to be merged (they are all CD34+) and a follwoing code will do the job;
```
GEX<-lapply(1:4, function(x){
dir=path[x]
tmp<-Read10X_h5(dir)
mat_raw <- tmp$`Gene Expression`
write_matrix_dir(mat = mat_raw, dir = dirname[x], overwrite = T, compress = F)
mat_bp<-open_matrix_dir(dir = dirname[x],buffer_size = 10000L)
CreateSeuratObject(counts = mat_bp, project = dirname[x])
  }
)
```
Comparing to the Seurat objs with nonBPCells matrix above, there are two additional lines that actually implements BPCells.  One is ```write_matrix_dir(mat = mat_raw, dir = dirname[x], overwrite = T, compress = F)``` which will create an indisk ```Formal class 'MatrixSubset' [package "BPCells"]```, followed by ```mat_bp<-open_matrix_dir(dir = dirname[x],buffer_size = 10000L)``` to access the indesk memory just created. Then all needed is to merge all the files and this could be done by ```merge`` command as decreived above.  

One immediate change in Seurat obj created by incorporating the BPcell matrix is the count matrix is no longer a ```Matrix obj``` but ```MatrixSubset obj``` as noted above, and normal Seurat downstream pipeline can be applied 

It is very simple to apply BPCells to the Seurat Objs with immediate benefits of having ```MatrixSubset objs```, most noticeably reduction in the file size ;
```
format(object.size(bp.merged.su), "GB")
[1] "0.3 Gb"
````
## Update on using BPCells
After playing with BPCells for a while, I decided to stop using BPCells for several reasons.  The most significant one that really bothered me was the lack of interoperability.  BPCells creates a marix, ```IterableMatrix object with class ColBindMatrices``` and this is not compatible  with other R and Python packages simply because they dont accept it.  Converting Seurat into SingleCellExperiment obj, for example, reuqires trasformation of ```IterableMatrix object with class ColBindMatrices``` into ```dgCMatrix```, and this poses some issues with current R, Rstuido and Positron (4.2.2), by giving an ***erros*** upon transposion of ```IterableMatrix object``` for some reasons that I have not figured out yet. Also, AnnData in Python will not accept this and I have to convert it to ```dgCMatrix``` or more preferably ```dgRMatrix```.  Instread of BPCells and having an headache and grave frustration with this, I decided to stick with ```Seurat's sketch```.   

