# Analysis of Combined CD34 and CD271 Samples 
## 1. Strategies for merging CD34 and CD271 Samples
Merging multiple ```Seurat files``` can be simply done by invoking ```merge``` command, however, this does not mean that any ```Seurat files``` can be merged easily.  For example, merging a seperate  ```CD34 Seurat file``` and ```CD271 Seurat file``` completed for preprocessing and analysis,  into a combined file does not work.  A process of untangling S3 structure of the multiple Seurat files with many slots gets too complicated to stich togher.  The best way to do this is to go from the beginning; simply importing various raw h5 files again and create a single Seruat object containing all the raw sample data. It could be an extremely large object but if BPCells or other means, the obj should be manageble in the local computer with sufficient memory for bioinformatics.  Or else, a big object like this may be better handled in the cloud.  

## 2. Importing 10X h5 files
Again the h5 file contains both GEX and ADT matrix, which can be imported by adding ```filename$`Gene Expression` or $`Antibody Capture` ``` repsectively. Following fuctions are exmaples of importing ```Gene Expression``` and ```Antibody Capture``` matrix; 
### Creating GEX Seurat Files;
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

### Creating ADT Seurat files;
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

## 3. SCTransform, PCA, Harmony Integration and UMAP Generation
These steps has become serioulsy routine, and in fact the latest Seurat V5 makes so much easier to run these.  
```
m.su.noMT<-SCTransform(m.su.noMT, assay = "RNA", new.assay.name = "SCT1", vars.to.regress = "percentMT", vst.flavor = "v2", verbose = T)
m.su.noMT<-RunPCA(m.su.noMT, assay = "SCT1", reduction.name = "SCT1.PCA", reduction.key = "SCT1PCA_", verbose = T, npcs = 50)
m.su.noMT<-RunHarmony(m.su.noMT, group.by.vars="orig.ident", reduction.use="SCT1.PCA", early_stop=F, plot_convergence=T, nclust=50, ncores=15, max_iter = 10, reduction.save = "SCT1.harmony")
m.su.noMT<-FindNeighbors(m.su.noMT, reduction = "SCT1.harmony", assay = "SCT1", graph.name = c("SCT_Harmony_knn", "SCT_Harmony_snn"))
m.su.noMT<-FindClusters(m.su.noMT, graph.name = "SCT_Harmony_snn", resolution = 6, algorithm = 4, cluster.name = "SCT_Harmony.snn.res6", method = "igraph")
```
I have to make a one comment on the last line begins with ```FindClusters```.  According to the paper, they identified **89** clusters which.  A argument/parameter of ```FindClusters``` that controls the cluster number is ```resolution```.  To get close to their number, I experimented with values of ```resulution```. 7 and 8 give over a hundred clusters while 6 yields 86.  I guess authors did not simply use a ```FindClusters``` function to derive the optimum number of clusters (86) but they also relied on the cell identifcation from the ADT assay to finally determine the cluster number. I will use **89 clusters** for now but trying to see how significant the differences are later. 
The next step, creating UMAP, is a bit tricky as the way I do.  Seurat could generate acceptable UMAP but quite often, shapes and distribution of points are not well separatable.  Seurate uses a R packages, ```uwot``` to generate UMAP with its default configuration (and Seurat's minor modifications). As with any R pacakges, ```uwot``` is customizable.   A real issue with custormerization, particularly with ```uwot``` is the fact that several parameters can control an ouput of the projection must be tweaked to generate an optimized projection but what is an the optimezed projection?  Interestingly, this has to be done by our own eyes.  Luckily, it is not necessaray to adjust every parameters but two parameters, ```a``` and ```b```.  I created a following function to test a range of values for ```a``` and ```b```.  
```
CreateSeuratUMAP<-function(suobj, reduction, a, b, neighbor, scale){
  tempfile(pattern = "model",tmpdir=tempdir(),fileext = ".uwot")->uwot_temp
  obj<-Embeddings(suobj,reduction = reduction)
  umap.obj<-umap2(obj, nn_method = "nndescent", n_threads = 19, a=a, b=b, fast_sgd = T, seed = 12345, ret_model = TRUE, n_neighbors = neighbor, n_components = 3, scale = scale)
  save_uwot(umap.obj, file = uwot_temp, unload = F)
  dt<-as.data.table(umap.obj$embedding)
  ggplot(dt,aes(V1, V2))+geom_point(alpha=0.5, shape=21)
  }
```
An example of the usage is as follows;
```
CreateSeuratUMAP(m.su.noMT, reduction = "SCT1.harmony", a = 0.9, b = 1, neighbor = 15, scale = "Z")
```
This function saves UMAP with particular values of ```a``` and ```b``` (```uwot_temp```) and draws a UMAP plot (by ggplot2).  Once the values of  ```a``` and ```b``` are found, they can be included in the Seurat obj as a ReducedDim obj, as follows;
```
m.su.noMT[["SCT1.UMAP"]]<-CreateDimReducObject(embeddings = model$embedding, assay = "SCT1", key = "SCT1UMAP_")
```


