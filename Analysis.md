## Data Preprocessing, Transformation, Dimension Reduction, and Annotation I
### Necessary steps to further process the Seurat obj for the downstrem analysis are written everywhere.  This section only covers those unique steps/processes that are unique to this data set, and are not well covered in the main stream Seurat tutorials.  

As always, the very first step of the standard preprocessing step is to filter out **low quality** cells, those with high MT transcript content and low RNA counts and features.  This step is well described everywhere and in my other analysis.  I will just mention briefly some of the key data preprocessing pointes regarding the data.  
### Creating  a ```SingleCellExperiment``` obj based on the ```BPCell count matrix``` for preprocessing
Unfortunately, this cant be done simple because the '''BPCell count matrix object``` is not compatible with how SingleCellExperiment stores the count (using ```Delayed Array```) matrix. However, there is a recently developed pakcage called ```BPCellsArray'''(https://github.com/Yunuuuu/BPCellsArray) makes it possible to utlize a ```BPCells``` like ***fuzzy in-memory operations*** in ```Delayed Array``` package format.  Like ```BPCells```, ```BPCellsArray``` provides the most of mathematical manipulations that are necessary to analye scRNAseq data in ```SingleCellsExperiment``` format.  
```
j<-JoinLayers(merge.su)
as(j@assays$RNA@layers$counts, "dgCMatrix") -> m
```
SinglCellExperimen still cant quite handle Seurat's layer concept yet, therefore, the Seurat obj must be joined first and extract the count matrix. 
Creating the ```BPCellArray``` count matirx is pretty simple;
```
path<-tempfile("BPCellsArray")
bitpacking_mat <- writeBPCellsDirMatrix(m, path = path)
bitpacking_mat
<36601 x 34742> sparse BPCellsMatrix object of type "double":
            AAACCCAAGATTCGCT-1_1 AAACCCAAGCGATTCT-1_1 AAACCCACAGCGACAA-1_1 ... TTTGTTGTCAGGACAG-1_4 TTTGTTGTCGGCACTG-1_4
MIR1302-2HG                    0                    0                    0   .                    0                    0
    FAM138A                    0                    0                    0   .                    0                    0
      OR4F5                    0                    0                    0   .                    0                    0
 AL627309.1                    0                    0                    0   .                    0                    0
 AL627309.3                    0                    0                    0   .                    0                    0
        ...                    .                    .                    .   .                    .                    .
 AC141272.1                    0                    0                    0   .                    0                    0
 AC023491.2                    0                    0                    0   .                    0                    0
 AC007325.1                    0                    0                    0   .                    0                    0
 AC007325.4                    0                    1                    1   .                    0                    0
 AC007325.2                    0                    0                    0   .                    0                    0

Storage Data type: double
Storage axis: col major

Queued Operations:
36601x34742 double, sparse: [seed] MatrixDir object
```
Astonishingly, size of ```bitpacking_mat``` is significantly smaller than the dgCMatrix, m.  
```
format(object.size(bitpacking_mat), "MB")
[1] "5.3 Mb"
format(object.size(m), "MB")
[1] "1925 Mb"
```
This diffrence also reflects the size differential ```SingleCellExperiment objs``` ;
```
sce1<-SingleCellExperiment(list(counts=m))
sce2<-SingleCellExperiment(list(counts=bitpacking_mat))
format(object.size(sce), "MB")
[1] "1932.2 Mb"
format(object.size(sce2), "MB")
[1] "12.5 Mb"
```
Interestingly, there is no significant outliers on the lower end;
```
isOutlier(sce2$nCount_RNA, type = "lower", nmads = 3, log = F) %>% summary()
   Mode   FALSE 
logical   34742
```
However there are some outlier cells at the higher end;
```
isOutlier(sce2$nCount_RNA, type = "higher", nmads = 3, log = F) %>% summary()
Mode   FALSE    TRUE 
logical   33312    1430
```
Let's look at some graphs representing some basic parameters;
![nCount](https://github.com/user-attachments/assets/2fc32daa-5a72-4bda-962e-96ed2c90e861)
Looks like data posted in GEO are preprocessed to some extend, already.  Other graphs of interest are;
![MTvsnCount](https://github.com/user-attachments/assets/bc5c5474-cb82-4ab2-8dd8-6dc3851b0ab1)
It appears that many cells with the high MT content (>75%) are not filtered out yet.  In the publication, cells with greater than 25% of MT content were reported to be removed.  
Also, a plot of ```nFeature_RNA``` vs ```nCount_RNA``` indicates that a majoiry of those ***RNA count outlier*** belong to the BM_27 samples.  
![Rplot](https://github.com/user-attachments/assets/947e0f31-05b8-4018-9a34-096dd797f20b)
Typically, the MT content is inversely correlated with RNA counts, which gives a good rationale for defining low quality cells  the high MT content and together with low RNA counts.  Removing cells with higher RNA features and counts are a bit tricky. A rationale for removing these cells shared by a dozens of publications is based on a vague assumption that these cells are most liekly doublelets. However this is not always the case; activated or proliferating cells naturally contain more RNA.  The best way to determine whether these cells are actually doubles is to run the doublet detection package, the best one in R being ```scDblFinder```.  

```scDblFinder``` detected potential doublets, and they appeart to be cells with higher RNA and Feature counts;
```
table(sce2$scDblFinder.class)
singlet doublet 
  30875    3867
```
![Rplot](https://github.com/user-attachments/assets/9075a075-a447-448a-b21f-cea3eb62b8b8)

Fortunately, there is no evidence of the presence of doublets disproportally distributed in the particular sample.  
```
table(sce2$scDblFinder.class, sce2$group)
         
           BF21  BM27  WF26  WM34
  singlet  6049  7210 10459  7157
  doublet   573   999  1239  1056
```
From purely the computtional stand point, these cells should be removed however, these may truely be singles with some biological significance.  It is wise not always rely completely on the mathmatical and staistical outputs for making a decision.  The best is to create a analysis pipeline with or without removing these doublets and see how they differ. 
But let's play with the singlet vs the doublet a bit longer before the next analysis step.  A follwing plot of the integrated UMAP from Seurat shows distribution  of the doublets.
![Rplot](https://github.com/user-attachments/assets/d1bb0b7b-3f36-41ab-8541-ee1537dc336a)
By looking at the cluster ditribution on the same UMAP, there are a few enriched doublets in the cluster 25 relative to other clusters.  
![Rplot01](https://github.com/user-attachments/assets/fe6d3e47-ae85-4389-888f-c4df1105fd6a)
This is more apparent by probing  a following ```scDblFinder.score``` plot. 
![Rplot02](https://github.com/user-attachments/assets/7d0fa774-9baa-44c6-9fd3-a47525264af4)
A dotpot of ```scDblFinder.score``` confirms this and also shows another cluster, 29 is also enriched with doublets. 
![Rplot03](https://github.com/user-attachments/assets/c81dd838-1193-4de9-8bc3-4b29738792ad)
Actual doublets counts among all the clusters indeed shows the cluster 25 has the most doublets.
```
    singlet doublet
  1     1314     219
  2     1331     192
  3     1077     217
  4     1175      50
  5     1101     120
  6     1130      88
  7     1159      46
  8     1080      78
  9     1071      56
  10    1039      84
  11    1082      26
  12     924     131
  13     969      75
  14     949      88
  15     951      58
  16     750     132
  17     783      91
  18     693     134
  19     721      84
  20     714      59
  21     669      72
  22     667      37
  23     612      53
  24     525     108
  25     264     369
  26     461     135
  27     409      75
  28     423      41
  29     266     146
  30     346      44
  31     354       5
  32     162      23
  33     127      27
  34      90      20
  35      74       2
```
These does not quite help deciding whether all the doublets should be removed or any doublets could be saved.  One factor that might augment the decisiion potentially is a differential gene expression of the doublets against singlets.  This might reveal the cellular nature of the doublets, for example, if a large portion of top DGEs are MT origin or Ribosome associated genes, these cells should well be removed.  








