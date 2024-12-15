## Analysis Step 1 -- Data Preprocessing, Transformation,  and Dimension Reduction
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
These does not quite help deciding whether all the doublets should be removed or any doublets could be saved.  One factor that might augment the decisiion potentially is a differential gene expression of the doublets against singlets.  This might reveal the cellular nature of the doublets, for example, if a large portion of top DGEs are MT origin or Ribosome associated genes, these cells should well be removed.  DGEs of doublets relative to singles could be easily done in ```Seurat```;
```
FindMarkers(j, ident.1 = doublet.name, assay = "RNA", logfc.threashold=0.5, test.use="wilcox", group.by = "group") %>% as.data.table(., keep.rowname=T) ->de.dt

de.dt
Index: <rn>
          rn         p_val avg_log2FC pct.1 pct.2     p_val_adj
      <char>         <num>      <num> <num> <num>         <num>
   1:   NAV1 2.803049e-175  0.2359848 0.486 0.220 1.025944e-170
   2:   PSD3 4.363296e-165 -0.1520466 0.544 0.272 1.597010e-160
   3: CDC25B 2.607588e-154  0.7004048 0.331 0.134 9.544034e-150
   4:  SPNS2 2.219335e-147  0.4199353 0.524 0.256 8.122989e-143
   5:   ASPM 2.538348e-146  0.6093897 0.461 0.228 9.290608e-142
  ---                                                          
8839: PTP4A2  8.662836e-01 -0.1119596 0.994 0.915  1.000000e+00
8840: NDFIP2  9.046675e-01 -1.5299037 0.018 0.018  1.000000e+00
8841: ATP5ME  9.213388e-01 -0.1384415 0.957 0.808  1.000000e+00
8842:   HBA2  9.709345e-01 -2.3858844 0.012 0.012  1.000000e+00
8843: RNASEK  9.772033e-01 -0.1158481 0.992 0.906  1.000000e+00
```
A total of almost 9000 DEGs !!  A volcano plot shows a majority of DGEs are low hanging fruits but top up and down regulated genes do not appear to be related to the apoptosis and no ribosome or mitochondria associated genes are found in this list of DEGs. Typically, low quality cells have hihg expressions of either or both ribosome and mitochondrial genes.  

![Rplot01](https://github.com/user-attachments/assets/ce793d67-d74d-4857-a310-95e2f1345228)

One last thing that might help assertain potential **phenotypes** is to perform functional profiling of this DEG set.  The quickest way to do this is to run a gGOSt analysis pipeline in the g:Profiler web site (https://biit.cs.ut.ee/gprofiler/gost).  I feed 250 DEGs ranked by p_val_adj (usually, the larger the better) and obtained a folloing output; 
![gProfiler_hsapiens_2024-11-21_01-25-31](https://github.com/user-attachments/assets/8dcd144c-83c5-4a28-b8cb-49a5c7f6613e)
 These compelling lines of analyses strongly suggest that these cells may not be junky doublets, though it is possible that those with extremely high RNA counts are indeed in the cycle and depending on where in the cycle, they may well be a doublet of progenitors.  At this point, these cells will be kept and included for the rest of downstream analysis.   Of note, Authors of this manuscirpt did not appear to run the doubelt discrimination.  
 
## Analysis Step 2 -- Annotation

Annotation is doen already by authors. In fact they really went distance with pain-staking steps to generate cell annotation though multiple statistical corrections including a ML step using XGboost, in order to nullify inconsistencies arised due to batch and integrative effects.  This whole process identified 89 distinct subsets in their combined samples, a total of 13 samples,  including titration controls.  Amount of the statistical scrunity is overwhelming and exhaustive. I am sure some of the reviwers asked them to run some of those statistical analysis, and there is no doubt about its completeness however, one bothersome question remaining for me is whether they went too far with this.  To publish a paper in Nature Med requires responses to extremely regorous questions from reviwers, and their extreme and rigrous staistical steps may be a testimony to that.  I can attest this from my experience dealing with a peer review process in Science. One reviewer had over 40 questions, which were extremely difficult to address.  Perhaps, this reviewer was purposely attempting to delay our acceptance  of our manuscript for his or her personal gain. Who knows.  
Gettig back to the cell annotation, I dont really have to run anything as mentioned. The cell annotation is provided by them, and all I have to do is to match UIDs in the cell annotation file to those of the transcriptome data in the Suerat format.  










