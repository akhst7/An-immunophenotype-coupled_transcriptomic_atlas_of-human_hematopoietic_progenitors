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
Let's look at UMAP created by this.  
![Rplot](https://github.com/user-attachments/assets/a3335dd0-51f6-490f-a172-9fa6a109dd48)

By looking at this overaly of all the sample onto the same UMAP space, ```Harmony``` worked really well.  A following plot shows the same integrated UMAP but with projection of cell annotation (86 distinct cell types).   
![Rplot01](https://github.com/user-attachments/assets/63f095ad-c968-46b3-8c89-e97c65e800f3)

Of note, there are 11992 cells with no annotations (NA in the legend). 

## RNA Velocity 
This is something I really wanted to try for this human HSC dataset.  In order to run the RNA Velocity analysis, counts of spliced and unspliced transcrpts must be generated by realighning fastqs from all the samples to the human genome.  There are several alingners such as Kallisto/BusTools and Salmon  that can do this but my go to aligner as of today is the ```Star``` aligner (https://github.com/alexdobin/STAR).  Typically, for 10x datasets this is very simple, and scRNAseq of this dataset is done by using 10x V3.1 kit which accroding to their protocal generates 28bp R1 (bar code + UMI) and 75 ~ 100bp R2 (cDNA).  However, the R1 and R2 fastq from this dataset contains 101bp.  Acoording to their method, they used PE150+10+10 with Illumina dual index for GEX sequencing. I am not sure what their rational to use this **unusal** approach as opposed to 10x recommended one, which is a bit confusing and the bottom line is that this makes the Star alignment a bit more difficult.  Besides, there are over 30 fastqs, each with 11GB at least and unless the HPC space is availableit, is not wise to download every relavant fastqs and run ```Star``` to get a Velocity count in a local computer.  
Anyhow, just take a quick look at two supposedly complementally R1 and R2 from the BF21_CD34+ sample;
```
(base) [~]$ seqkit stat /Volumes/MySlateDrive/GSM7836870/SRR26369695_1.fastq                        
file                                                  format  type     num_seqs         sum_len  min_len  avg_len  max_len
/Volumes/MySlateDrive/GSM7836870/SRR26369695_1.fastq  FASTQ   DNA   139,892,272  14,129,119,472      101      101      101
(base) [~]$ seqkit stat /Volumes/MySlateDrive/GSM7836870/SRR26369695_2.fastq
file                                                  format  type     num_seqs         sum_len  min_len  avg_len  max_len
/Volumes/MySlateDrive/GSM7836870/SRR26369695_2.fastq  FASTQ   DNA   139,892,272  14,129,119,472      101      101      101
```
R1 and R2 are identical in length, meaning they may well be complementary afterall they were generated by PE150 but they are not.  First few sequneces of R1 look as follows;
```
head /Volumes/MySlateDrive/GSM7836870/SRR26369695_1.fastq                                                                                   
@SRR26369695.1 A00587:951:HWLG3DSX5:1:1101:1434:1000 length=101
CTTTACCAGACAACATTCCGGGCCCCCAAAACACCAGCCCAACCGCACACATATTTTTAATCCCGGCCCCCACGGTTCCCCAGCAGGAGGTCCGGGGTTTG
+SRR26369695.1 A00587:951:HWLG3DSX5:1:1101:1434:1000 length=101
,FFFFFFFFFFFFFFFFFFFFFFFFFFF,:F:F,F,:,FF,,:,F:F,,,::,,F:FF,,:,:F,:F,FFF,FFFF,,:::,F:,FF,:::,F,F,:,:,:
@SRR26369695.2 A00587:951:HWLG3DSX5:1:1101:2682:1000 length=101
TTGCGTCAGTCAGCCCTTAGGCATTAATTATATACCTACAGAGAAGCGCAGCGGTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+SRR26369695.2 A00587:951:HWLG3DSX5:1:1101:2682:1000 length=101
FFFFFFFFFFFFFFFFFFFFFFFFFFFF,:,F,F,::,::,,,F,,F,,,F,,F,,:,FFF:,FFFF:FF:F:F,F,FF:F,FF,FF,:FFFFF,,:,,,F
@SRR26369695.3 A00587:951:HWLG3DSX5:1:1101:3007:1000 length=101
GTCGCGAGTTGTCATGCGCTGACGTAGCGGATGATCACAGCGCAAATAAGTACCCCACCACACAAAATACACAGACAGGTGTTTCAATAATCACTCTCACA
```
And R2;
```
(base) [~]$ head /Volumes/MySlateDrive/GSM7836870/SRR26369695_2.fastq                                                                                   
@SRR26369695.1 A00587:951:HWLG3DSX5:1:1101:1434:1000 length=101
GGTGGATCACTGGAAGCCAGGAGTTTGAGACAAGCCTGGCCAACATGGTGAAACCCCGTCTCTACTAAAAATAACAAAAATTAGACGGGCACGGTGGTGTG
+SRR26369695.1 A00587:951:HWLG3DSX5:1:1101:1434:1000 length=101
FFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:F,FFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@SRR26369695.2 A00587:951:HWLG3DSX5:1:1101:2682:1000 length=101
GCCCCCATTCGTATAATAATTACATCACAAGACGTCTTGCACTCATGAGCTGTCCCCACATTAGGCTTAAAAACAGATGCAATTCCCGGACGTCTAAACCA
+SRR26369695.2 A00587:951:HWLG3DSX5:1:1101:2682:1000 length=101
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFF:FFFFFFFFFFFFFFFFFFF
@SRR26369695.3 A00587:951:HWLG3DSX5:1:1101:3007:1000 length=101
TGGAATCAACACCACCGAGCTCTGTGGGAAAAAAGAAAAACCTGCTCCCTTCGCTCTGCTGGAAGCTGGAGGGTGCTAGGCCCCTGTGTAGTAGTGCATAG
```
Reverse complement of R2 should look exactly like R1 or vice versa but this is not the case.  
```
(base) [~]$ seqkit head -n 1 /Volumes/MySlateDrive/GSM7836870/SRR26369695_2.fastq |seqkit seq -r -p -t dna                                                                   
[INFO] when flag -t (--seq-type) given, flag -v (--validate-seq) is automatically switched on
@SRR26369695.1 A00587:951:HWLG3DSX5:1:1101:1434:1000 length=101
CACACCACCGTGCCCGTCTAATTTTTGTTATTTTTAGTAGAGACGGGGTTTCACCATGTTGGCCAGGCTTGTCTCAAACTCCTGGCTTCCAGTGATCCACC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFF,F:FFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFF
```
This R2 rev-comp  does not look complemtary to R1.  Further examination shows that R1 contains the bar codes;
![Screenshot 2025-01-03 at 9 24 56 PM](https://github.com/user-attachments/assets/caf8a572-4f30-419c-9dbe-a9c6847594b5)
```BF21CD34bc.txt``` is a list of entire barcodes used in the BF21_CD34+ sample. (Of note, the first sequence of R1 does not contain any barcodes !!)
According to 10x, the set of next 12 sequences following the bar code corresponds to UMI,and supposedly if the sequencing is done PE, the following set of sequences corresponds to cDNA. However they are not complementary to R2 whatsoever.  
![Screenshot 2025-01-03 at 9 36 13 PM](https://github.com/user-attachments/assets/de379f5a-915f-4eee-936a-1e4bcbaa2231)
Do reverse comped R2 sequences contain any BF21_CD34+ barcodes ? Answer is no.  
![Screenshot 2025-01-03 at 9 50 51 PM](https://github.com/user-attachments/assets/45eb189c-483a-467a-821c-47079859e281)

I really hate to go after RNA Velocity, at this point.  I will find out more about this and update this but the most important goal at this point is to create a SingleCellExperiment Object based on the current Seurat object for developing a human hematopoietic cell annotation referece for scRNAseq.  

### Update 
I posted a question to the first author, Xuan Zhang  of this article and he suprisingly got back to me promptly.  According to Xuan, the library was created according to a 10x V3.1 protocol, and fastqs were generated by CellRanger, meaing that R1 and R2 are prototypical fastqs that contian 28bps (16bp barcode + 12bp UMI) and 100bp cDNA repectively.  Then, R1 needs to be truncated to retain the first 28 bps. Also a strand infomation of both R1 and R2 needs to be corrected as well.  A following is an example of how to clean up R1 and R2 fastqs prior to running STARsolo to obtain spliced and unspliced counts. Of note, there are so many different ways to accomplish the same goal.  It goes all the way back from the begining, downloading fastqs from GEO.  

```
A script to clean up fastqs for 10x V3 compatibility 
prefetch GSM7836870 -O /Volumes/Temp/BF21_CD34 #downloadig necessary SRA corresponding a BF21 CD34 sample. 
fd . -t f /Volumes/MySlateDrive/BF21_CD34 >sra.file 
parallel -j 20 --bar  -vv 'fasterq-dump -e 20 -m 50G -v  -p -t /Volumes/Bioinformatics -O /Volumes/MyslateDrive/BF21_CD34 {}' :::: sra.file
fd . -d 1  -t f   /Volumes/MySlateDrive/BF21_CD34 |rg _1 >R1.list                                
fd . -d 1  -t f   /Volumes/MySlateDrive/BF21_CD34 |rg _2 >R2.list
parallel -j 20 -vv 'gsed -i "s/^\+.*/+/g" {}' ::: *.fastq
parallel -j 20 -vv 'gsed -i "2~4 s/\(^.\{28\}\)\(.\{73\}\)/\1/" {}' ::: *_1.fastq
parallel -j 20 -vv 'gsed -i "4~4 s/\(^.\{28\}\)\(.\{73\}\)/\1/" {}' ::: *_1.fastq
parallel -j 20 -vv 'gsed -i "s/length=101/length=28/g" {}' ::: *_1.fastq
parallel -j 20 -vv --bar 'pigz --fast -p20 -v {}' ::: *.* # STARsolo prefers gz files but I may be wrong about this one. 
```
Once a correct version of R1 and R2 are created, STARsolo can be run as follows; 
```
#!/opt/homebrew/bin/zsh

index=/Volumes/MySlateDrive/star_genome_index/index_human_Gencode_GRCh38_p14_47
whitelist=/Volumes/Bioinformatics/star_genome_index/3M-february-2018.txt
r1=/Volumes/Bioinformatics/temp/SRR26369695_1.fastq.gz
r2=/Volumes/Bioinformatics/temp/SRR26369695_2.fastq.gz
r3=/Volumes/Bioinformatics/temp/SRR26369696_1.fastq.gz
r4=/Volumes/Bioinformatics/temp/SRR26369696_2.fastq.gz
r5=/Volumes/Bioinformatics/temp/SRR26369697_1.fastq.gz
r6=/Volumes/Bioinformatics/temp/SRR26369697_2.fastq.gz
r7=/Volumes/Bioinformatics/temp/SRR26369698_1.fastq.gz
r8=/Volumes/Bioinformatics/temp/SRR26369698_2.fastq.gz

STAR --genomeDir $index \
--runThreadN 50 \
--readFilesCommand 'pigz -c -d' \
--readFilesIn $r2,$r4,$6,$8 $r1,$r3,$5,$7 \
--soloCBwhitelist $whitelist \
--soloType CB_UMI_Simple \
--soloCBstart 1 \
--soloCBlen 16 \
--soloUMIstart 17 \
--soloUMIlen 12 \
--soloStrand Forward \
--twopassMode Basic \
--clipAdapterType CellRanger4 \
--soloMultiMappers EM \
--soloFeatures Gene GeneFull Velocyto \
--outFilterMultimapNmax 1 \
--outFilterScoreMin 30 \
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
--soloUMIfiltering MultiGeneUMI_CR \
--soloUMIdedup 1MM_CR \
--soloCellFilter  EmptyDrops_CR \
--soloCellReadStats Standard \
--outReadsUnmapped Fastx \
--outFileNamePrefix /Volumes/Bioinformatics/BF21_34 \
--outSAMattributes NH HI CR CB CY UR UY sM \
--outSAMtype BAM SortedByCoordinate
```
STARSolo creates many ouputs but the most important and relevant files are in a ```xxxSolo.out/Velocyte/Filtered``` folder.  
![Screenshot 2025-01-21 at 3 03 59 PM](https://github.com/user-attachments/assets/548d0ed9-8d8b-495c-bed3-d9faea40f090)

All the files necessary to create either a ```Seurat``` or ```Scanpy/Scvelo```,  including count, barcode and cell files are in this folder.   Notably, showen here is an example using fastqs from the one sample,  BF21_CD34.  For entirey, fastqs from all the samples are needed to be realigned to get spliced vs unspecified counts, and this requires enormous computing power and resources which even with my MacStudio it will be difficult, although my MacStudio with M1Ultra takes only 20min to complete STARSolo run, suggesting it would only take a day to finish realigning all the samples.  My issues is not so much of the lack of computing power but data storage.  My system is short of the memory space enough to host all the Solo ouputs.  

One might argue that rather than realigning fastqs from the individual sample, which requires multiple Solo runs,  concatenating all the fastqs and run Solo in a single run. However, an issue of doing this is that some of the same barcodes were used repeatedly in several samples thereby, fastqs with the same bar code yet from different samples will be counted as the same.  The only way this woud work is if @RG has the SM (group) which can be used to distinguish samples in BAM files.  

Actual sequence of the RNA velocity analysis will be covered in another more dedicated git page for the RNA velometry and psuedotime analysis. 

## Establishing the scRNAseq Cell Annotation Ref for the Human Hematopoietic Tissue
This dataset comes with cell annotations based on the ADT(CITE-Seq) and transcriptomic data, then why further cell annotation step is necessary? It turns out that only **11992** cells out of 67889 (in this example buy a lot more in the actual published data) are annotated.  These unannotated cells could be annotated by using annotated 67889 cell as a training set for annotating the test set, 11992 unannotated cells.  This actually gives a great opportunity to test how well an annotation algorithm performs.  

There are many R as well as Python packages that can create farily precise and accurate cell annotations.  More recent ones take advantage of ML and AI (e.g. ChatGPT) algorithm for feature selection, a feature in this case is cell annotation/type, while more established ones typically use the simpler approach and they work just as well as ones with ML and AI.  I use a R package called, ```SingleR```.  Since the original version, (https://www.nature.com/articles/s41590-018-0276-y) was published back in 2019,  several improvements by a current maintainer, Aaron Lun have been added.  ```SingleR``` depends on ```SummerizedExperiment```, thus needs ```SingleCellExperiment``` format, and therefore, the count matrix in the ```Seurat``` has to be imported to the ```SingleCellExperiment``` obj.  Since ```SingleCellExperiment``` does not accept the ```IterableMatrix``` count matrix in the ```Seurat``` obj, it has to be converted first to ```dgCMatrix```.  Also, ```SingleCellExperiment``` does not accept **layered data** as in ```Seurat V5``` obj, and these layered matrix must be joined to form a single count matrix.  
```
joined<-JoinLayers(m.su.noMT)
```
Then, the joined ```IterableMatrix count matrix``` will be converted to ```dgCMatrix```;
```
as(johined@layers$counts, "dgCMatrix")-> m
```
There is one another step needed prior to assembling a ```SingleCellExperiment``` obj, that is to transfer meta data from the ```Seurat``` to ```SingleCellExperiment``` obj.  ```SingleCellExperiment``` does not accept a typiecal dataframe but ```DataFrame``` by ```S4Vector```.  
```
DF<-DataFrame(m.su.noMT[[]])
sce<-SingleCellExperiment(list(counts=m), colData=DF)

> sce
class: SingleCellExperiment 
dim: 36601 67889 
metadata(0):
assays(1): counts
rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
rowData names(0):
colnames(67889): AAACCCAAGGCCCACT-1_1 AAACCCACAAGAGCTG-1_1 ... TTTGTTGTCAGGACAG-1_8 TTTGTTGTCGGCACTG-1_8
colData names(13): orig.ident nCount_RNA ... Level3M Level3R
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):

> colData(sce)
DataFrame with 67889 rows and 13 columns
                             orig.ident nCount_RNA nFeature_RNA percentMT percentRB nCount_SCT1 nFeature_SCT1 SCT_Harmony.snn.res8 seurat_clusters
                            <character>  <numeric>    <numeric> <numeric> <numeric>   <numeric>     <integer>             <factor>        <factor>
AAACCCAAGGCCCACT-1_1 BPCells/BF21-CD271       4877         2288   4.44946   13.1843       11280          2679                   24              37
AAACCCACAAGAGCTG-1_1 BPCells/BF21-CD271      11472         3720   3.52162   22.1409       12095          3717                   71              68
AAACCCACAATACCTG-1_1 BPCells/BF21-CD271       8030         2599   6.77460   30.8095       11988          2605                   2               2 
AAACCCACAATGACCT-1_1 BPCells/BF21-CD271       4657         1866   5.06764   22.8688       11566          2451                   53              61
AAACCCACATATGAAG-1_1 BPCells/BF21-CD271       7923         2629   2.90294   33.6994       11964          2644                   32              22
...                                 ...        ...          ...       ...       ...         ...           ...                  ...             ...
TTTGTTGAGATAGCAT-1_8  BPCells/WM34-CD34      25079         5037   3.44511   38.4425       26305          5036                  3                7 
TTTGTTGAGTCTCGTA-1_8  BPCells/WM34-CD34      51039         7424   3.79514   26.8324       27906          6564                  100              3 
TTTGTTGCATTATGCG-1_8  BPCells/WM34-CD34      31558         5500   1.49249   38.0918       27824          5500                  86               14
TTTGTTGTCAGGACAG-1_8  BPCells/WM34-CD34      25834         4978   2.57413   36.8700       26468          4977                  86               14
TTTGTTGTCGGCACTG-1_8  BPCells/WM34-CD34       4409         1827   2.44954   19.4829       25060          4735                  35               48
                     SCT_Harmony.snn.res7 SCT_Harmony.snn.res6        Level3M      Level3R
                                 <factor>             <factor>    <character>  <character>
AAACCCAAGGCCCACT-1_1                   31                   37         BMCP-1         BMCP
AAACCCACAAGAGCTG-1_1                   69                   68             NA           NA
AAACCCACAATACCTG-1_1                   11                   2              NA           NA
AAACCCACAATGACCT-1_1                   55                   61             NA           NA
AAACCCACATATGAAG-1_1                   33                   22             NA           NA
...                                   ...                  ...            ...          ...
TTTGTTGAGATAGCAT-1_8                   2                    7           HSC-1        HSC-1
TTTGTTGAGTCTCGTA-1_8                   94                   3  MultiLin-GMP-2 MultiLin-GMP
TTTGTTGCATTATGCG-1_8                   60                   14       immNeu-1       immNeu
TTTGTTGTCAGGACAG-1_8                   60                   14             NA           NA
TTTGTTGTCGGCACTG-1_8                   80                   48        Pro-B-2      Pro-B-2
```
One of the improvements of ```SingleR``` implemented by LT allows ```SingleR``` to run parallelley by utilizing a ```BiocParallel``` package.  To take advantage of this to speed things up, an parallel instruction in the form of ```BiocParallel``` object needs to be created;
```
library(BiocParallel)
MulticoreParam(18, RNGseed = 1234, progressbar = T)-> bp

> bp
class: MulticoreParam
  bpisup: FALSE; bpnworkers: 18; bptasks: 2147483647; bpjobname: BPJOB
  bplog: FALSE; bpthreshold: INFO; bpstopOnError: TRUE
  bpRNGseed: 1234; bptimeout: NA; bpprogressbar: TRUE
  bpexportglobals: TRUE; bpexportvariables: FALSE; bpforceGC: FALSE
  bpfallback: TRUE
  bplogdir: NA
  bpresultdir: NA
  cluster type: FORK
```
The last step, prior to running ```SingleR```, is remove ```NA``` from a ```Level3M``` ref column of the colData.  
```
sce.new.Nona<-sce[, !is.na(sce$Level3M)]

> sce.new.noNa
class: SingleCellExperiment 
dim: 36601 55897 
metadata(0):
assays(2): counts logcounts
rownames(36601): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
rowData names(0):
colnames(55897): AAACCCAAGGCCCACT-1_1 AAACCCAGTAACGCGA-1_1 ... TTTGTTGCATTATGCG-1_8 TTTGTTGTCGGCACTG-1_8
colData names(14): orig.ident nCount_RNA ... Level3R sizeFactor
reducedDimNames(0):
mainExpName: NULL
altExpNames(0):

> is.na(sce.new.noNa$Level3M) %>% sum()
[1] 0
```
Now, ```sce.new.Nona``` can be used as a cell annotation ref in the following ```SingleR``` line;
```
test.anno<-SingleR(test = sce.na, ref = sce.new.noNa, labels = sce.new.noNa$Level3M, de.method = "wilcox", BPPARAM = bp)
```





