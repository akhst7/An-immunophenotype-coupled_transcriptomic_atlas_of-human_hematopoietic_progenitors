# Analysis of Combined CD34 and CD271 Samples 
## 1. Strategies for merging CD34 and CD271 Samples
Merging multiple ```Seurat files``` can be simply done by invoking ```merge``` command, however, this does not mean that any ```Seurat files``` can be merged easily.  For example, merging a seperate  ```CD34 Seurat file``` and ```CD271 Seurat file``` completed for preprocessing and analysis,  into a combined file does not work.  A process of untangling S3 structure of the multiple Seurat files with many slots gets too complicated to stich togher.  The best way to do this is to go from the beginning; simply importing various raw h5 files again and create a single Seruat object containing all the raw sample data. It could be an extremely large object but if BPCells or other means, the obj should be manageble in the local computer with sufficient memory for bioinformatics.  Or else, a big object like this may be better handled in the cloud.  

## 2. Importing 10X h5 files
Again the h5 file contains both GEX and ADT matrix, which can be imported by adding ```filename$`Gene Expression` or $`Antibody Capture` ``` repsectively. Following fuctions are exmaples of importing ```Gene Expression``` and ```Antibody Capture``` matrix; 
# Creating GEX Seurat Files;
```
fs::dir_ls("../An immunophenotype-coupled_transcriptomic_atlas_of human_hematopoietic_progenitors/GSE245108/") %>% .[grep("CD34|CD271",.)] ->path #GSE245108 is a folder containing all the h5 files
cellanno[!grep("TNC", Group), .(unique(Group) %>% str_replace(., "-", "_"))] -> dirname
fs::dir_create(paste0("BPCells/", dirname)) #These are directories that ```BPCells``` stores on disk fuzzy data. 
CreateSeuratBPCellMatrix<- function(x){
  f=path[x]
  tmp<-Read10X_h5(f)
  mat_raw <- tmp$`Gene Expression`
  write_matrix_dir(mat = mat_raw, dir = dirname[x], overwrite = T, compress = F)
  mat_bp<-open_matrix_dir(dir = dirname[x],buffer_size = 10000L)
  CreateSeuratObject(counts = mat_bp, project = dirname[x])
}
l<-lapply(1:8, CreateSeuratBPCellMatrix)
names(l)<-dirname
merged.su<-merge(l$`BF21_CD271`, l[2:8]) #This will merge 8 Seurate objects that are created with the BPCell matrix.  In Seurat V5, all the individual Seurat objs are stored as layers. 
```

# Creating ADT Seurat files;
```
adt<-lapply(1:8, function(x){
  N<<-1+N
  dir=path[x]
  tmp<-Read10X_h5(dir)
  a<-tmp$`Antibody Capture`
  colnames(a)<-paste0(colnames(a),"_",N)
  CreateAssay5Object(a)
}
)
names(adt)<-dirname
adt.m<-merge(adt$`BF21-CD271`, adt[2:8], add.cell.ids=c("_1", "_2", "_3", "_4","_5", "_6", "_7", "_8")) #add.cell.ids will add "_x" to the barcode as a suffix. Interestingly, this is not necessarily in merging "GEX" matrix)
```
As mentioned in the previous section, the ADT matrix are considerable smaller than GEX matrix.  Rows only consistes of 142 which corresponds to the number of antibodies in the cocktail used, and there is no need to call ```BPCells```.  Also, the ADT matrix will be incorporated as a ```Seurat Assay5 Object``` rather than the ```Seurat Object```. Once this is done, ADT V5Assay objs can be incorporated into GEX Seurat Obj as an assay.  
```
merged.su[["ADT"]]<-m.adt
```
One thing to lookout for is the number of columns (cell number) between merged.adt and merged.su.  If they dont agree, Seurat will complain. 

# SCTransform, Harmony Integration and UMAP Generation
