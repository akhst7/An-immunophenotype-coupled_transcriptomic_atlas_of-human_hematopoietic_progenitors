## Data Preprocessing, Transformation, Dimension Reduction, and Annotation I
### Necessary steps to further process the Seurat obj for the downstrem analysis are written everywhere.  This section only covers those unique steps/processes that are unique to this data set, and are not well covered in the main stream Seurat tutorials.  

The very first step of preprocessing step is to filter out **low quality** cells, those with high MT transcript content and low RNA counts and features.  It appears that 
