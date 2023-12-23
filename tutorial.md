# `FlowCT`: A semi-automated workflow for deconvolution of immunophenotypic data and objective reporting on large datasets
#### _(by: Cirino Botta, Catarina Maia and Juan-Jose Garcés)_

`FlowCT` should be used as an analysis pipeline to retrieve all available information from _FCS_ files. It is scalable, being able to analyze from thousands to millions of cells and from tenth to hundreds of markers (depending on computational power availability). It includes different algorithms for batch removal and for cell clustering. Additionally, by using the `SingleCellExperiment` object class, it has wide compatibility with existing packages designed for single-cell RNAseq analysis. 

Here we propose an ideal pipeline and we will explain, step by step the most important functions of the package by applying it to a full example. This tutorial is based on a serie of paired samples (peripheral blood [PB] and bone marrow [BM] from the same subjects) from a total of 10 patients (data available [here](https://flowrepository.org/id/RvFrgzHozvhHBC1p5vfTppnIDqiln8bv3QZqgX5OCcyLN6Kipdd5OM5eIun6IMNK)). For a deeper explanation and additional options in each function please carefully read package documentation (`?function_name`).

> ***Before starting***: All the _FCS_ files should include compensation before being analyzed with `FlowCT`. This step should be done with commonly available flow cytometry software. Additionally, all the samples should have been stained with the same panel. It is not important that the name of the marker is included within the _FCS_ nor the channels have the same order.

> Our pipeline accepts both unprocessed and pre-processed (i.e., manual doublet removal or population pre-selection) _FCS_ files.

## 1. Metadata and `SingleCellObject` (`SCE`) preparation
Choose the working directory which include the _FCSs_ to be analyzed: `setwd("your/path/with/fcs")`

Verify that all _FCSs_ has the same marker names. As already stated, this is an important step. We included a function to perform a quick check and fix (providing the same markers order) of eventual inconsistencies...


```{r}
library(FlowCT)
unify.FCSheaders(directory = ".", pattern = "fcs", fix = F)
```

```
#              names     	 freq	   length
# 1 FSC-A:NA, ...truncated..., V450-A:CD27, V500-A:CD45RA               27        13
# 2 FSC-A:NA, ...truncated..., Pacific Blue-A:CD27, AmCyan-A:CD45RA      3        13
```


...if name inconsistencies are found, the option `fix` should be changed to `TRUE` and select the desired order within those showed by the function (if not provided, the most frequent one will be choosen). All the files with inconsistencies will be replaced by fixed ones (with the _fixed_ preffix) while a copy of the original files will be maintained in the _./original_files/_ folder (within the working directory). If some files have a different number of markers, they will be discarded in a new folder called _./discarded_files_because_diffs/_.

### Metadata matrix construction

The metadata matrix includes the most important information regarding the files to be analyzed. We have two possibilities to generate this matrix:

1. Outside R, by preparing a file which include all the information of the files to be analyzed (with Excel or text editors); 

2. inside R, by manually creating the matrix. 

The latter one is the easiest approach if there is need/possibility of quickly adding or removing files from the analysis at different time-points. To do that, we chose to use a specific code for the name of each _FCS_, by including all the relevant information separated by a specific character. We use the code _condition_patientID_other.fcs_ that let us to calculate a new metadata matrix very quickly every time we add (or remove) new files to (from) our folder. 

```{r}
(filenames <- list.files(pattern = "fcs", path = "."))

md <- data.frame(filename = filenames, #mandatory
     sample_id = 1:length(filenames),
     condition = sapply(filenames,function(x) strsplit(x, split = "_|\\.")[[1]][1]),
     patient_id = sapply(filenames,function(x) strsplit(x, split = "_|\\.")[[1]][2]))

# note: it could be useful to transform some variables to factor
md$condition <- as.factor(md$condition)
```

### `SCE` object generation


The `SCE` class is a lightweight Bioconductor container initially born for storing and manipulating single-cell genomics data (see [here](http://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) a more detailed description). However, we found extremely useful its use with flow cytometry data. To do that, original fluorescence data and transformed data were used in the place of RNA counts and logcounts.

We designed a function (`fcs.SCE`) which include several options to optimize the generation of the SCE from a bunch of _FCS_ files.

```{r}
fcs <- fcs.SCE(directory = ".", pattern = "fcs", metadata = md, events = 10000, transf.cofactor = 500, project.name = "ImmTh")
```


- Taking into account that this pipeline has been firstly designed for the analysis of large number of files with an high events number from clinical trials, a <u>downsampling</u> (reduction of total events per file to be considered in the analysis) step is usually performed (`events` option). This step is optional and strongly dependent on the computational power of the system used for downstream analyses. We here selected 10000 events for each file.

- Furthermore, we added a <u>transformation</u> option to deal with the bad representation of negative data values with digital cytometer: the Log transformation from the cytometer leads to compression of data against the axes, poor visual representation of low intensity or unstained populations and errors in subsequent automatic gating/clustering algorithms. Accordingly, we used `arcsinh` transformation by using a cofactor of 500 (which should be optimized according to data in the `trans.cofactor` option; e.g. _FCSs_ from the Aurora spectral cytometer would need a cofactor of 4000). Both the transformed and the non-transformed matrix are stored in the `SCE` object within the `"transformed"`and `"raw"` assays, respectively. 


### Marker names replacement and identification of group of markers for further analyses


In this step we will proceed to rename markers in a format more suitable for representation and to avoid conflicting symbols in R (e.g., "FITC-A: CD62L" will be transformed in "CD62L"). Additionally, we will select groups of markers to be used in downstream analyses (please note that even if not correct, we included FSC-A and SSC-A within the group of surface markers. This has been done to simplify the use of these markers in the clustering procedure). 

```{r}
marker.names(fcs) #see current channel names
new_names <- c("FSC_A","FSC_H","SSC_A","SSC_H", "CD62L", "CXCR3", "CD8", "CCR4", 
      "CCR6","CD4", "CD27", "CD45RA")
fcs <- marker.names(fcs, new.names = new_names)

# markers selection
surface_markers <- c("FSC_A", "SSC_A","CD62L", "CXCR3", "CD8", "CCR4", "CCR6",
      "CD4", "CD27", "CD45RA")
physical_markers <- c("FSC_H", "SSC_H")
```

### Automatic quality control and doublet removal


Once we have a `SCE` built on events-reduced/names-checked/compensated files, we can start the quality check. To do that, we designed a one-step function which will perform an automatic detection/removal of unwanted anomalies and doublets. The algorithm removes suspected anomalies deriving from changes in the flow rate, signal acquisition instability and outliers/margin events in the dynamic range (based on the [`flowAI`](https://www.bioconductor.org/packages/release/bioc/html/flowAI.html) package). Next, the function eliminates the doublets by using a code that automatically design a density-based ellipsoid gate in the FSC-A/FSC-H dotplot and use the 95th percentile as limit for singlet cells (code has been  written by [Droumeva R.](https://github.com/LennonLab/flowcytometry/blob/master/bin/support_functions.R)).

```{r}
fcs <- qc.and.removeDoublets(fcs, physical.markers = c("FSC_A", "FSC_H", "SSC_A", 
      "SSC_H"), return.fcs = F)
# with return.fcs = T, high QC FCS will be saved in a separate folder
```

### Detection and removal of batch-effect


Despite the great efforts in the standardization of flow cytometry protocols, biological and technical variability can be reduced but not completely abrogated. To deal with this problem we offer four different approaches which follow two different philosophies. 

- The first two, `gaussNorm` and `warp` methods (included within the [`flowStats`](https://www.bioconductor.org/packages/release/bioc/html/flowStats.html) package), intent to reduce the technical variation by "shifting" and aligning the area of high local density (landmarks) in each channel. This is a per <u>channel modification</u> and could be applied with the function `marker.normalization`.

```{r}
fcs <- marker.normalization(fcs.SCE = fcs, marker = c("CCR6", "CCR4"), method = "warp")
```

- The second approach, more holistic, relies on the possibility of "aligning" clusters of cells and modify the expression matrix accordingly to the new position. This wrapper function (`batch.removal`) includes the single-cell (RNAseq) methods `Seurat` and `harmony` (please, refer to the original development sites, [[1]](https://satijalab.org/seurat/v3.0/integration.html) and [[2]](https://github.com/immunogenomics/harmony]), for additional information). 

```{r}
fcs <- batch.removal(fcs, method = "harmony", batch = "patient_id")
```


> While `gaussNorm` and `warp` methods can easily handle millions of events, the computational power required for both `harmony` and (especially) `Seurat` usually limits the number of events that could be processed to hundreds of thousands. 

From an operative point of view, the first step to perform in this process is to look at the distribution of peaks in each channel and identify the parameter(s) to be aligned.

```{r}
multidensity(fcs.SCE = fcs, assay.i = "transformed", subsampling = 100)
```


![Figure 1. Multidensity's output](https://i.ibb.co/nzXJ3kk/Slide1.png)

This function present two additional options that could be useful in certain scenarios:

- It is possible to limit the parameters to be shown through the option `show.markers = c("param1", "param2")`. 

- It is possible to select if histograms should be shown as "stacked" or "overlapped" histogram. The latter simply rely on the decision of a cut-off (`ridgeline.lim = x`): if your _FCSs_ are beyond this value they will be represented as "stacked" otherwise we will see an "overlay" histogram.

By looking at the histograms we identified many markers which are not perfectly "aligned". Anyway, by taking into account that about 50% of the whole cellularity is composed by granulocytes which could interfere with the batch removal procedure, we decided to apply this step lately to the lymphocytes population only. 

## 2. Descriptive and exploratory analysis 
### Total cell count per group


At this point, it could be useful to check for possible confounding variables that may affect further downstream analysis, such as important differences in the total number of cells between conditions. It is possible to choose the variables to use to group patients (`x.axis` option). 

```{r}
cell.count.bx(fcs.SCE = fcs, x.axis = "condition", assay.i = "transformed")
```


As reported in the Figure 2A, we have a significant difference in the cell number of PB and BM samples, which could be expected due to the better quality of blood samples... but it will barely affect our subsequent analysis.

### Principal component analysis (PCA) and heatmap on median marker expression per sample


Dimensionality reduction (DR) algorithms aim to maintain the structure of the data (i.e., the distance between the different subpopulations or patients in our case) while reducing the dimensions needed to observe differences (by jointing similar objects and separating different ones). PCA is one of these techniques and perform DR by projecting data to new coordinates preserving the variance. 

Mimicking a transcriptomic approach, [Nowicka M _et al._](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5473464/) proposed to use the median marker expression calculated over all cells to highlight similarities/dissimilarities between samples. We implemented this approach with `median.dr` and `median.heatmap` functions. We also added the possibility to use `kmeans` to cluster patients according to their characteristics (ellipses will be eventually generated, in the DR, and identified clusters will be available for downstream analyses). 

```{r}
cell.count.bx(fcs.SCE = fcs, x.axis = "condition")
median.dr(fcs, color.by = "filename", label.by = "condition", assay.i = "transformed", 
            num.k = NULL) # number of kmeans to calculate (and drawing)

median.heatmap(fcs, assay.i = "transformed", 
      not.metadata = c("sample_id", "filename")) # hide some features from annotation
```

Of note, both representations well discriminate PB from BM samples. This could be due to the different cellular composition that could affect median fluorescence values (e.g., PB does not include erythroid precursors which are negatives for all markers). This point will be clarified by performing a quick look at "single-cell level". 

![Figure 2. Initial exploratory analysis](https://i.ibb.co/X4D5Gc9/Slide2.png)

### PCA and heatmap on single cell expression


We decided to perform DR (both PCA and UMAP [Uniform Manifold Approximation and Projection]) and to generate a clustered heatmap at single-cell level (by subsampling 500 cells from each patient). 

```{r}
fcs500 <- sub.samples(fcs.SCE = fcs, subsampling = 500)
fcs500 <- dim.reduction(fcs500, dr.method = c("UMAP", "pca"), assay.i = "transformed")

dr.plotting(fcs500, plot.dr = "PCA", color.by = "condition", size = 0.8, 
      assay.i = "transformed", colors = c("#c03728", "#919c4c")) 
dr.plotting(fcs500, plot.dr = "UMAP", color.by = "condition", size = 0.8, 
      assay.i = "transformed", colors = c("#c03728", "#919c4c")) 

sc.heatmap(fcs500, markers.to.use = surface_markers, assay.i = "transformed")
```

![Figure 3. Single-cell functions' outputs](https://i.ibb.co/qsWHmrp/Slide3.png)



As reported in Figure 3, single-cell PCA (A) demonstrated an almost superimposable distribution between samples irrespective of the grouping according to tissue origin. Even better, the UMAP (B) reports an almost identical distribution of cells derived from PB or BM, with the exception of a specific cluster of cells more abundant in the BM, presenting a complete negativity for the expression of all surface markers (as reported in the heatmap, C)... as we mentioned, they are likely to be erythroblasts, physiologically absent in PB. 

## 3. Cell Clustering
As already stated before, analyzing datasets of flow cytometry data with 8 or more colors (for exploratory/discovery purpose) with traditional gating strategies is absolutely not easy, time-consuming and difficult to reproduce. Additionally, our capability to work with two-dimensional plots (or rarely with three-dimensional ones) only, renders quite impossible to keep an overview of all markers behavior on selected cell types, increasing the risk of missing valuable information. 

To overcome these limits, several methods of unsupervised cell clustering have been developed. Among the different available methods, we included in our pipeline 4 different clustering approaches: `FlowSOM`, `Phenograph`, `Seurat` and the recently introduced `PARC`. Each of these methods have its drawbacks that should be carefully evaluated at the time of clustering method decision. It is often a good idea to compare clustering results before choosing with which one proceed to next steps. `FlowSOM` is the quickest one, while `Phenograph` and `Seurat` are the most precise. `PARC` is the newest one and it seems to balance precision and speed. 

The function to perform cell clustering is `clustering.flow` and include different options that could be added to each specific method (`parc`, `som`, `Seurat`, `phenograph`). Here we will use `FlowSOM` to quickly identify and select lymphocytes for subsequent, more detailed, analyses. Within the `num.k` option it is possible to select the number of final clusters we expect to identify (please take into account that the `FlowSOM` algorithm works better if you "overcluster" your dataset (i.e., you select a number of clusters higher than the populations you expect to identify).

`fcs <- clustering.flow(fcs, method = "SOM", assay.i = "transformed", num.k =15)`

### Dimensionality reduction

The application of dimensionality reduction algorithms (already described previously) to single-cell expression data leads to the production of beautiful easy-to-interpret images that give us the idea of distance between groups of cells. Additionally, each cell could be colored according to the relative expression of each marker analyzed or to the previously identified cluster. Due to the high computational power necessary to perform such analyses, it is often convenient to reduce the cells number to be analyzed for each sample just before running the dimensionality reduction code. To be noted that by performing this step at this point, we do not affect the clustering process that has already been performed on the whole dataset. 

For our example, taking into account that our aim is to identify blood and BM macro-populations, we decided to perform dimensionality reduction (PCA and UMAP) on 500 cells per sample. By side, our `dim.reduction` function offer the possibility to calculate the t-stochastic neighbour embedding (t-SNE), currently less used due to its limitations such as the high computational demand, long computational times and difficulty to effectively reduce large datasets or a very large number of dimensions (usually in single-cell -omics, where dimensions for each cell are thousands, t-SNE is calculated on PCA-reduced dimensions to reduce the computational load). The obtained reduced dimensions could be plotted in a `ggplot2` graph by using the function `dr.plotting as` showed below.

> In order to generate a palette of colors, a possible option is (check [another color schemas](https://www.datanovia.com/en/wp-content/uploads/dn-tutorials/ggplot2/figures/0101-rcolorbrewer-palette-rcolorbrewer-palettes-colorblind-friendly-1.png)): 
> `color_clusters <- RColorBrewer::brewer.pal(n = 12, name = "Set1")`

```{r}
fcs500 <- sub.samples(fcs.SCE = fcs, subsampling = 500)
fcs500 <- dim.reduction(fcs500, dr.method = c("UMAP", "pca"), assay.i = "transformed")

dr.plotting(fcs500, plot.dr = "UMAP", assay.i = "transformed", colors = color_clusters, 
				color.by = "SOM.k15") # UMAP according to FlowSOM clusters

# color.by is not expecified, all marker will be shown (faceted)
dr.plotting(fcs500, plot.dr = "UMAP", assay.i = "transformed")

median.heatmap(fcs.SCE = fcs, assay.i = "transformed", 
               cell.clusters = "SOM.k15") # heatmap grouped according to FlowSOM clusters
```

![](https://i.ibb.co/SmdPKcc/Slide4.png "Figure 4")


> ***Tips***: the `dr.plotting` function allows for figure faceting according to a selected variable. This could be done by adding the option `facet.by` to the function. 


## 4. Clusters identification, _FCS_ exporting and clusters merging
Once that all the processes of clustering and dimensional reduction have been completed, we can proceed to the next step which include the identification of each cluster as a concrete population. This could be done by exporting an "enriched" _FCS_ and analyzing it outside of _R_, or within _R_ by combining UMAP, heatmaps and dot-plots images.

### _FCS_ exporting and clusters analysis

We can export our down-sampled dataset in a new _FCS_ file containing for each cell: the DR coordinates (e.g., UMAP1, UMAP2, PCA1, PCA2); the number of clusters according to selected method; other metadata variables and the fluorescence intensity value from the <u>untrasformed</u> matrix. This file could be loaded in any flow cytometry software and subsequently analyzed to identify the population to which each cluster belongs. At this point the analysis becomes familiar to researchers with experience in flow cytometry, indeed, what should be done now is to identify for each cluster, by manual gating, the cell population which it belongs. 

![](https://i.ibb.co/BBxrLNV/Slide5.png "Figure 5. Manual annotation on external software")

In Figure 5B, we drawn a gate around the cluster number 1 and by looking at the distribution in bidimensional dot plots (physical parameters and CD8-SSC dot plots are reported as the most useful) we could conclude that these cells belong to the "granulocytes" group. In our dataset, we assigned to 21 out of 40 clusters the name "lymphocytes", and in Figure 5C are reported the reliability of these results. Specifically, by manual gating we were able to identify the 22.76% of total cellularity as lymphocytes, a result almost identical to the one obtained by merging the 10 clusters identified by `FlowSOM` that we recognized as lymphocytes (23.56%). By fusing the XYZ lymphocytes clusters, we observed that the algorithm includes erroneously in these groups about a 3% of cells (0.71% of total cellularity) that are not lymphocytes; on the other side, no CD4+ or CD8+ lymphocytes could be found in the remaining clusters. This error, that could already be considered acceptable, could be eventually reduced in further subclustering analyses focused on lymphocytes group only (see after).   

### Cluster identification within R

Within `FlowCT`, we included a specific function, `flowplot`, that allows users to depict flow cytometry dot-plot and histograms to help to identify the cell population included within the cluster. To be noted that it is possible to select the clusters to be plotted, which additionally improves the possibility to correctly identify each cluster. 

```{r}
flowplot(fcs500, x.axis = c("CD62L", "CD4", "CD4"), y.axis = c("CD45RA", "CD8", 
      "SSC_H"), color.by = "SOM.k25", assay.i = "transformed")
```

> Of note, it is always better to plot expression data from the original, "untouched", matrix as compared to the batch-removed one.  

![](https://i.ibb.co/v355h33/Slide6.png "Figure 6. Flowplot's output")


### Clusters merging

Once all clusters have been correctly identified, it is possible to produce final images reporting the merged clusters as reduced derived clusters with the chosen name. In our dataset the 25 clusters were divided into the subsequent macro-populations: 18/25 lymphocytes, 3/25 granulocytes, 2/25 erythroblasts/debris, 1/25 eosinophils, 1/25 monocytes. It is possible to upload the new cluster names in 2 ways: by entering them within the script, or by providing a spreadsheet including old and new names (we will use the first approach here and the second subsequently).

```{r}
# prepare data.frame with pairs: cluster-population (it could be done on the fly)
replacedata <- data.frame(original_cluster = levels(fcs$SOM.k25))
replacedata$new_cluster <- c(rep("granulocytes", 3), "lymphocytes", "eosinophils", 
      "granulocytes", rep("lymphocytes", 3), "erythroblasts/debris", "monocytes", 
      "lymphocytes", "lymphocytes", "erythroblasts/debris", "monocytes", 
      rep("granulocytes", 10))

fcs500$SOM_named <- clusters.rename(fcs500$SOM.k25, cluster = replacedata$original_cluster, 
      name = replacedata$new_cluster)
fcs$SOM_named <- clusters.rename(fcs$SOM.k25, cluster = replacedata$original_cluster, 
      name = replacedata$new_cluster)

dr.plotting(fcs500, plot.dr = "UMAP", color.by = "SOM_named", assay.i = "transformed", 
      colors = color_clusters)
median.heatmap(fcs.SCE = fcs500, assay.i = "transformed", cell.clusters = "SOM_named", 
      colors = color_clusters)
```

![](https://i.ibb.co/X2nXJNh/Slide7.png "Figure 7. Clusters merging")

### Subclustering

At this point, the next step is to perform a deeper analysis on specific subpopulation identified with the `FlowSOM` clustering followed by cluster merging. Our panel have been designed to evaluate T cells population, and, accordingly, we decided to perform a further analysis on lymphocytes. In this view, it could be useful (especially in term of speed), in specific settings, to prepare a new SCE object including all the events belonging to lymphocytes population.

```{r}
fcsL <- fcs[,fcs$SOM_named == "lymphocytes"]
```

Again, we calculated the cells number for each condition (to be noted how, obviously, PB samples present a higher number of lymphocytes as compared to BM ones) and calculate an heatmap on the median expression for each marker.

```{r}
cell.count.bx(fcsL, assay.i = "transformed", x.axis = "condition")
median.heatmap(fcsL, assay.i = "transformed", not.metadata = c("sample_id", "filename"))
```

![](https://i.ibb.co/f2G3LDW/Slide8.png "Figure 8. Subsampling a cell population")

### Batch removal

At this point, we decided to proceed with the `Seurat` method for batch removal. It is important to disclose, within the function, the variable that include the batch (in our case all patients have been acquired in different days and this should be considered the most important batch; we assume that BM and PB of the same patients, being stained and acquired simultaneously by the same person are not affected by any batch) and the name of the new assay that will be generated. Additionally, further options specific for each method could be passed to the function by using the argument `seurat.params` or `harmony.params`.

> ***Optional***: the `Harmony`'s method is often a good and quick alternative to Seurat, however, in our experience, we observed Seurat as a better method when handling less than 12 markers.
> ***Important***: this could be a looooooooooong step.

```{r}
fcsL <- batch.removal(fcsL, method = "seurat", batch = "patient_id", 
      new.matrix.name = "normalized", seurat.params = list(reduction = "cca", dims = 1:10))
```

We can now look at the histograms before and after batch correction:

```{r}
multidensity(fcs.SCE = fcsL, assay.i = "transformed", subsampling = 100)
multidensity(fcs.SCE = fcsL, assay.i = "normalized", subsampling = 100)
```

![](https://i.ibb.co/f24bxsS/Slide9.png "Figure 9. Compare normalized data")

We can further perform a check in a dimensionality reduced plot to observe the presence and eventual disappearance of batch effect after batch.removal function.

```{r}
fcsL <- dim.reduction(fcsL, dr.method = c("PCA", "UMAP"),assay.i = "normalized")
fcsL_old <- dim.reduction(fcsL, dr.method = c("PCA", "UMAP"),assay.i = "transformed")

dr.plotting(fcsL_old, plot.dr = "UMAP", assay.i = "transformed", 
   colors = color_clusters, size = 0.3,color.by = "patient_id")
dr.plotting(fcsL, plot.dr = "UMAP", assay.i = "normalized", colors = color_clusters, 
   size = 0.3, color.by = "patient_id")
```

![](https://i.ibb.co/PC1vTdH/Slide10.png "Figure 10. Compare normalized data (UMAP)")

> ***Tips***: Due to the low number of dimensions used, in some analysis `Seurat` could produce infinite or NAs values within the expression matrix which could impairs some subsequent functions. If this happens a workaround could be to assign a "0" value to `NA`s (this should be carefully took into account) and  a "-99999" or "99999" to positive or negative infinite values.
> 
> ```{r}
> assay(fcsL, "normalized")[is.infinite(assay(fcsL, "normalized"))] <- -99999
> assay(fcsL, "normalized")[is.na(assay(fcsL, "normalized"))] <- 0
> ```

## 5. Clustering on lymphocytes only
### Cell clustering and dimensional reduction

As done before, the next steps will be to apply a similar workflow to the new reduced object to obtain detailed information about all the lymphocyte subpopulations present in our dataset. We will then perform a clustering with the `PARC` algorithm and will apply it to the previously calculated PCA and UMAP. Lastly, we will perform a cluster identification (by either exporting a new FCS containing only the lymphocytes subpopulation or by taking advantages of FlowCT tools) to correctly characterize each subpopulation. 

```{r}
fcsL <- clustering.flow(fcsL, method = "PARC", assay.i = "normalized")

dr.plotting(fcsL, plot.dr = "UMAP",color.by = "PARC", assay.i = "normalized", 
   colors = color_clusters, size = 0.3)
median.heatmap(fcs.SCE = fcsL, assay.i = "normalized", cell.clusters = fcsL$PARC)

dr.plotting(fcsL, plot.dr = "UMAP", assay.i = "normalized", colors = color_clusters, 
   size = 0.3)

export.metaFCS(fcs.SCE = fcsL, output.name = "subclust1")
```

![](https://i.ibb.co/vwyDYLn/Slide11.png "Figure 11. Subclustering analysis")

### Identification and representation of different lymphocytes subpopulations

Finally, in this last step of the clustering process, we could assign the population name to each cluster and update the clustering representations to show us the correct names. As shown in Figure 11, we were able to identify 28 different lymphocytes subpopulations (represented in UMAP and heatmap on median marker expression). All subpopulations are well separated and recognizable. 

```{r}
replacedataL <- readxl::read_excel("Th_new_clustering.xlsx", col_types = "text")
fcsL$PARC_L_named <- clusters.rename(fcsL$PARC, cluster = replacedataL$original_cluster, 
   name = replacedataL$new_cluster)

dr.plotting(fcsL, plot.dr = "UMAP", color.by = "PARC_L_named", colors = color_clusters, 
   size = 0.3)
median.heatmap(fcs.SCE = fcsL, assay.i = "normalized", cell.clusters = fcsL$PARC_L_named)

dr.plotting(fcsL, plot.dr = "UMAP", color.by = "PARC_L_named", colors = color_clusters, 
   size = 0.3, facet.by = "condition")
```

![](https://i.ibb.co/m4Q79nm/Slide12.png "Figure 12. Subclustering cell asignment")

Another useful tool for check marker's implication on each cell population: 

```{r}
## differential dotplots
dotplot.DE(fcs.SCE = fcsL_rm, markers.to.use = surface_markers, 
   clusters.named = "PARC_L_named")
```

![](https://i.ibb.co/Q9qLcSL/Slide14.png "Figure 13. Differential dotplot")

## 6. Tools for statistical analysis and representation of the complexity
### Phylogenic tree

To have a better idea of distribution and phenotype of each population we represented our results in a phylogenic tree taking advantage of the recently developed [`ggtree`](https://www.bioconductor.org/packages/release/bioc/html/ggtree.html) package. We designed a wrapper function, `circ.tree`, which initially let us show the nodes to select as population dividers (see below) and then let us color them accordingly. The heatmap associated to each cell population reports the median expression of each marker, while the dimension of the circle is proportional to the relative abundance of the population itself: Naïve CD4 cells are the most abundant T population in our dataset. 

```{r}
circ.tree(fcs.SCE = fcsL, cell.clusters = "PARC_L_named", nodes = "display")
circ.tree(fcs.SCE = fcsL, cell.clusters = "PARC_L_named", nodes = c(29, 36,31))
```

![](https://i.ibb.co/6wTfDty/Slide12.png "Figure 14. Circular clustering")

### Populations abundance visualization and export of table of results 

The most useful thing that we can obtain from this analysis, in a clinical trial perspective, is a table including all the abundance (in term of both absolute count and percentage) or median fluorescence of each identified population for each patient. Such table indeed, could be subsequently used in any statistical software (even out of _R_) and population proportions could be correlated with clinical parameters such as time to endpoint, response to treatment, etc.  

Within `FlowCT`, the relative abundance of each population could be reported in a barplot (including a bar dedicated to each sample and grouped for tissue of origin), and as a table, including all results, readily available for downstream use.

```{r}
prop_tableL <- barplot.cell.pops(fcs.SCE = fcsL_rm, cell.clusters = "PARC_L_named", 
   count.by = "sample_id", facet.by = "condition", return.mode = "percentage")
dataset <- data.frame(md, as.data.frame.matrix(t(prop_tableL)))
write.table(dataset, file = "results.txt", sep = "\t")


med <- median.values(fcs.SCE = fcsL_rm) # median values (from non-transformed data)
dataset2 <- data.frame(md, med)
write.table(dataset2, file = "results_median.txt", sep = "\t")
```

![](https://i.ibb.co/PtLkmKc/Slide13.png "Figure 15. Cell percetajes distribution")


> ***Tips***: Often, clusters which include debris or unknown populations are identified. It is important (and useful) to remove these clusters for statistical analysis. Please, be aware of what you are removing.
>
> ```{r}
> delete_pops <- c("Unc", "debris")
> fcsL_rm <- remove.pop(fcsL, clusters.named = "PARC_L_named", population = delete_pops)
> ```


### Evaluation of correlation and comparison analyses  

Subsequently, we performed a correlation analysis between PB and BM to investigate potential differences in the distribution of each subpopulation. To this aim, we provide a function, `corplot.condition`, which allow us to decide which condition to correlate and which parameters should be investigated for correlation. Interestingly, we found a direct correlation for almost all populations between PB and BM.

```{r}
corplot.conditions(fcs.SCE = fcsL_rm, cell.clusters = "PARC_L_named", 
   condition = "condition")
```

![](https://i.ibb.co/X30tB9r/Slide15.png "Figure 16. Correlation between two conditions")

Lastly, we compared the abundance of each population according to tissue of origin. We provide different methods to visualize comparison results (and respective significance levels) which are available with the following code.

```{r}
## boxplot
boxplot.cell.clustering(fcs.SCE = fcsL_rm, cell.clusters = "PARC_L_named", 
   facet = T, return.stats = F, plot.only.sig = c(T, 0.05), facet.free.scale = "free")
sig <- boxplot.cell.clustering(fcs.SCE = fcsL_rm, cell.clusters = "PARC_L_named", 
   return.stats = T)

## Dumbbell plot
dumbPlot(fcs.SCE = fcsL_rm, cell.clusters = "PARC_L_named", return.stats = F, 
   condition.column = "condition")

## parallel plot
parallel.plot(fcs.SCE = fcsL_rm, cell.clusters = "PARC_L_named", return.stats = F, 
   condition.column = "condition")
```

![](https://i.ibb.co/b35KQCb/Slide16.png "Figure 17. Differential plots")


### Machine learning for immune populations selection

Between all machine learning (ML) modeling advantages, feature selection is extremelly useful for identify those populations with a strong implication with our event of interest (for example, disease progression). Thus, we also parepared a (very basic) wrapper including some ML approaches for "playing" with our identified immune populations: 

- [biosigner](https://bioconductor.org/packages/release/bioc/html/biosigner.html), which helps to assess the relevance of variables through different <u>binary</u> (i.e., they don't take into account the survival times).
- [randomForestSRC](https://kogalur.github.io/randomForestSRC/theory.html)... the classical random forest approach, also being able to use survival times for regression (and not only for classfication, like the previous one).
- [SurvBoost](https://github.com/EmilyLMorris/survBoost) helps to select more relevant variables (i.e., immune populations) through gradient boosting. Note: for installing this package, please, go to indicated link (it's not available within CRAN repository).

They use either percentages (or absolute number of events) or categorized values (according different cutoffs calculations through `FlowCT::pop.cutoff`) depending of `cutoff.type`'s values. Because each ML method has a lot of internal parameters for tunning them, you can also include them within our function, always as list. Let's see some minimal examples with `SurvBoost`:

```{r}
# add survival data
s <- readxl::read_excel("survival_data.xlsx")
colData(fcsL_rm) <- merge(colData(fcsL_rm), s, by = "patient_id")
head(colData(fcsL_rm))
```

```
# DataFrame with 6 rows and 9 columns
#    patient_id            filename sample_id condition     cell_id  PARC       PARC_L_named       PFS  censored
#   <character>         <character> <integer>  <factor> <character>  <factor>   <factor>      <numeric> <numeric>
# 1           1 BM_1_017294.red.fcs         1        BM         1.4  12   CD8_TM_CXCR3p_CD27p        10         1
# 2           1 BM_1_017294.red.fcs         1        BM         1.9  12   CD8_TM_CXCR3p_CD27p        10         1
# 3           1 BM_1_017294.red.fcs         1        BM        1.17  6    CD4_N                      10         1
# 4           1 BM_1_017294.red.fcs         1        BM        1.30  5    B_cells                    10         1
# 5           1 BM_1_017294.red.fcs         1        BM        1.32  1    NK_cells                   10         1
# 6           1 BM_1_017294.red.fcs         1        BM        1.33  10   CD8_EM_CXCR3p_CD27p        10         1
```

```{r}
ml1 <- prog.pop.selection(fcs.SCE = fcs, cell.clusters = "PARC_L_named", 
                          time.var = "PFS", event.var = "censored", 
                          cutoff.type = "terciles",
                          # these are internal parameters for controlling survboost
                          method = "survboost", 
                          method.params = list(rate = 0.1, 
                            control_method = "num_selected", 
                            control_parameter = list(num_select = 5)))

ml1bis <- prog.pop.selection(fcs.SCE = fcs, cell.clusters = "PARC_L_named", 
                          time.var = "PFS", event.var = "censored", 
                          cutoff.type = "none", # percentages are used as they are
                          method = "survboost", 
                          method.params = list(rate = 0.1, 
                            control_method = "likelihood"))
```

![](https://i.ibb.co/KFdbcC2/Slide23-ML.png "Figure 18. ML approaches")

And, of course, you can also use the classical `predict` function when using training and validation datasets (in the case of `SurvBoost` is slightly different, just for knowing):

```{r}
# select the 70% of your samples for training (set.seed, IMPORTANT!!!)
set.seed(333)
train_idx <- sample(unique(fcs$filename), length(unique(fcs$filename))*0.7)

ml2 <- prog.pop.selection(fcs.SCE = fcs, cell.clusters = "PARC_L_named", 
                          time.var = "PFS", event.var = "censored", train.index = train_idx, 
                          cutoff.type = "maxstat", method = "survboost", 
                          method.params = list(rate = 0.1, 
                            control_method = "num_selected", 
                            control_parameter = list(num_select = 3)), 
                          return.ML.object = T) #TRUE, for later predict

SurvBoost::predict.boosting(ml2$ML.object, newdata = ml2$survival.data$data$test)
```


## Appendix A: interaction with other packages for single cell analysis (`scater` or `Seurat`)
We designed `FlowCT` keeping in mind the possibility to interconnect it with pipelines already available for single cell analysis. Within the _R_ environment, `scater` and `Seurat` represent 2 of the most widely used toolkits and each of them includes tools that could be useful for flow cytometry analysis. By using the `SCE` structure for storing flow cytometry data, we made transition within the different packages easy.

### Scater

To use `scater` commands with `FlowCT` objects, the easiest thing to do is to create a new assay named `"logcounts"` and transfer there the ones we would like to use: 

```{r}
library(scater) 

assay(fcsL, "logcounts")<-assay(fcsL, "transformed")
```

We can now use all the useful utilities for data visualization included in `scater` such as: `plotExpression`, which could be used in combination with `flowplot` to correctly assign cluster identities; the dimensionality reduction functions derived from `plotReducedDim` including `plotUMAP`;  `ggcells` function, which perfectly interact with `ggplot2` for custom graph such as density plot:

```{r}
plotExpression(fcsL, rownames(fcsL)[9], x = "condition", colour_by="condition")+ 
	scale_fill_manual(values = color_clusters) + theme(legend.position="none")
plotExpression(fcsL, rownames(fcsL)[7], x = "PARC", colour_by="PARC") + 
	scale_fill_manual(values = color_clusters) + theme(legend.position="none")

plotUMAP(fcsL, colour_by = "CD8")

ggcells(fcsL, mapping=aes(x=UMAP.1, y=UMAP.2, colour=CD4)) + 
	geom_point() + stat_density_2d() + 
	facet_wrap(~condition) + 
	scale_colour_distiller(direction = 1) + theme_minimal()
```

![](https://i.ibb.co/jz26n3G/Slide17.png "Figure 19. Scater integration")

Furthermore, if needed, it is possible to calculate dimensionality reductions from `scater` commands including `runUMAP` or `runPCA`.

### Seurat

The conversion from a `SCE` experiment to a `Seurat` object has been recently made "painless". Indeed, the `Seurat` package provides a function dedicated to this kind of "migration"; importantly, it should be indicated the assays that should be passed to the counts and data slots:

```{r}
library(seurat)

fcsL_Se <- as.Seurat(fcsL, counts = "raw", data = "transformed")
fcsL_Se <- ScaleData(fcsL_Se, features = rownames(fcsL_Se))
```

At this point, all the functions commonly used for single cell RNAseq analysis are available for flow cytometry data. Of note, `Seurat` batch removal and clustering functions are already included in `FlowCT`, but it is possible to perform also these analysis by using the `Seurat` object. All dimensional reductions could be calculated and plotted within `Seurat` and visualizing tools including heatmaps could be easily implemented in flow cytometry pipeline. 

We can use the UMAP representation to visualize cell populations:

```{r}
DimPlot(fcsL_Se, reduction = "UMAP", group.by = "PARC_L_named", repel = T, label = T) +
	scale_color_manual(values = color_clusters) + theme(legend.position="none")
```

Or visualize marker expression, with the possibility to set upper and lower cut-off to improve visualization or to "blend" two markers to visualize double expressing cells

```{r}
FeaturePlot(fcsL_Se, features = c("CCR4"), pt.size=0.1, min.cutoff = "q10", 
  split.by = "condition", max.cutoff = "q90") +
  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "PuBu")[2:9])

FeaturePlot(fcsL_Se, features = c("CD4", "CXCR3"), blend = TRUE, blend.threshold = 0.1)
```

![](https://i.ibb.co/ZNg8T2k/Slide18.png "Figure 20. Seurat integration")

Furthermore, we can identify the markers which are significantly over-represented in each population and plot them in a heatmap:

```{r}
Idents(object = fcsL_Se) <- "PARC_L_named"
fcsL.markers <- FindAllMarkers(fcsL_Se, only.pos = TRUE, min.pct = 0.25, 
   logfc.threshold = 0.25)
top5 <- fcsL.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

DoHeatmap(subset(fcsL_Se, downsample = 100), features = top5$gene, size = 3) + 
   NoLegend()
```

![](https://i.ibb.co/c2JcSpW/Slide19.png "Figure 21. Seurat integration (heatmap)")

Or we can use all others visualization tools present in seurat such as dotplot and ridgeplot

```{r}
DotPlot(fcsL_Se, features = rownames(fcsL_Se), cols = c("blue", "red"), 
   dot.scale = 8, split.by = "condition") + RotatedAxis()

RidgePlot(fcsL_Se, features = c("CD8","CD4"), ncol = 2)+ NoLegend()
```

![](https://i.ibb.co/k1wt5ZF/Slide21.png "Figure 22. Seurat integration (misc viz)")

## Appendix B: pseudotime analysis (`Slingshot`/`Monocle`)
The objective of pseudotime analysis is to measure the relative progression of each cell along a specific process (such as maturation or differentiation) allowing us to understand the (pseudo)temporal behaviour without explicit time series. Of note, this analysis is possible when elements within the dataset behave asynchronously and each is at a different stage of "evolution". The algorithm, by creating a relative ordering of the cell, can define a series of "states" that constitute a trajectory for the process of interest. 

To explain this step, we will use a dataset including the CD8 population from the 24 colors example. Batch effect has been removed by using the canonical correlation algorithm from `Seurat` package and the corrected expression matrix has been saved under the `"Seurat"` assay. As first step, we will calculate a `DiffusionMap` by using the `scater` package:

```{r}
fcsL <- runDiffusionMap(fcsL, exprs_values = "Seurat")

plotDiffusionMap(fcsL, colour_by = 'PARC') +  scale_fill_manual(values  = color_clusters)

# if marker expression on DiffusionMap is planned to be visualized, it is important, 
# as previously reported, to create the "logcounts" assay within the SCE object
assay(fcsL, "logcounts")<-assay(fcsL, "transformed")
```

We then selected the clusters including the CD8 populations only:

```{r}
fcsL_pseudo <- fcsL[, !fcsL$PARC=="5" & !fcsL$PARC=="0" & !fcsL$PARC=="8" & 
   !fcsL$PARC=="12" & !fcsL$PARC=="14"& !fcsL$PARC=="13"& !fcsL$PARC=="17"& 
   !fcsL$PARC=="15"& !fcsL$PARC=="23"& !fcsL$PARC=="20"& !fcsL$PARC=="29"& 
   !fcsL$PARC=="27"& !fcsL$PARC=="16"]

plotDiffusionMap(fcsL_pseudo,colour_by = 'PARC') + 
   scale_fill_manual(values = color_clusters)
```

To calculate the "pseudotime" trajectory we take advantage of the `slingshot` package. As reported by the authors, the purpose of this package "is to use clusters of cells to uncover global structure and convert this structure into smooth lineages represented by one-dimensional variables, called "pseudotime". It is important to note that slingshot does not require previous knowledge on cell cluster evolution (i.e., it does not need a "start" and "end" cluster). Due to computational power necessary for subsequent steps, we subsampled our dataset just to take 3000 CD8 lymphocytes from each sample. Next, we performed pseudotime calculation on `PARC` clusters and `DiffusionMap` DR by using the `slingshot` command.

```{r}
library(slingshot)
library(RColorBrewer) 

fcsL_sling <- sub.samples(fcs.SCE = fcsL_pseudo, subsampling = 3000)
plotDiffusionMap(fcsL_pseudo, colour_by = 'PARC') + 
   scale_fill_manual(values = color_clusters)

sim <- slingshot(fcsL_sling, clusterLabels = 'PARC', reducedDim = 'DiffusionMap')
```

Pseudotime calculation will be stored within the `sim` object. It could be easily moved within the original object for downstream analysis. 

```{r}
fcsL_sling$Pseudotime<-sim$slingPseudotime_1
plotDiffusionMap(fcsL_sling[,!is.na(fcsL_sling$Pseudotime)], colour_by="Pseudotime")
```

Lastly, we can take advantage of the package monocle 2 to depict a "smoothed" heatmap showing all the significant changes that happens in markers involved in pseudotime progression. To do that, we will convert our `SCE` object in a `Seurat` object and then in a `monocle` object (a working conversion tool from `SCE` to `monocle` is not available at this moment) 

```{r}
library(monocle)
set.seed(333)

fcsL_sling_SEU <- as.Seurat(fcsL_sling, counts = "raw", data = "transformed")
fcsL_sling_mon <- as.CellDataSet(subset(fcsL_sling_SEU, downsample = 150), assay = "RNA")
```

Next, we will create a new object that only contains the cells for which a pseudotime value has been calculated. This will be then used for the identification of genes that are significantly modulated along "progression". Those genes will be then represented in an heatmap with the `plot_pseudotime_heatmap` function. Of note, it is possible to automatically group markers that are modulated in the same way to help understanding the biological substrate of pseudotime. Detailed information on `slingshot` and `monocle` are available on the developer's websites. 

```{r}
cds_subset<-fcsL_sling_mon[,!is.na(fcsL_sling_mon$Pseudotime)]

diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.05))

plot_pseudotime_heatmap(cds_subset[sig_gene_names[c(1:5,7:14)],], num_clusters = 3, 
   cores = 2, show_rownames = T, return_heatmap = TRUE, norm_method = "vstExprs")
```

![](https://i.ibb.co/cY0y4R3/Slide22.png "Figure 23. Pseudotime application")
