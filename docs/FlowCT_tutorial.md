## FlowCT (R package)
### A semi-automated workflow for deconvolution of immunophenotypic data and objective reporting on large datasets

Hello! Thanks for interesting you in FlowCT!

This is a very generic tutorial showing some possibilites of this R package. For any additional information or check another analysis appoximations, you can see our paper in [xxxxx](wwww.web.com) and the help for each function of this package.

# Adjust this part with general pipeline below!
# Explain data used for tutorial (PB/BM)
So let's go! Pipeline is structurated in multiple parts, some internals and other externals:
1. Preprocess of FCS files, ie. to unify (florofore) channels names and set them in the same order. Normaly, it depends of computation power available, user also wants to reduce the number of events to avoid computer crashing or to facilitate compensation of multiple FCS.
2. External compensation of FCS files (if needed) to correct overlap between channels. Whatever flow cytometry software can be used, in our case we use [Infinicyt<sup>TM</sup>](http://www.infinicyt.com/).
3. Import data from all FCS compensated files, check uniformity of marker names (eg. CCR4 and not CD194), quality control, doublet removal and transformation/normalization.
4. Descriptive analysis through heatmaps, boxplots, PCA... at general and single-cell level.
5. Cell clustering using SOM (Self-Organizing Map), PCA, tSNE and UMAP and respective visualizations.
6. New global FCS generation and manually identification of generated clusters in the previous step. Later, assign and merge of identified cell populations.
7. (Optional) Subclustering and repeat the analysis on a desired specific population.
8. Exporting table for statistical analysis, either inside R or another statistic software.

# Add references and URLs for used packages (?)
#### 0) Prepare FCS files
Try to avoid whitespace character in the filenames of your FCS files and rename them with an informative name (see below the reason, *spoiler*: metadata generation is based on filenames).
We've written a simple code to perform this task (on [Windows](https://github.com/jgarces02/FlowCT/tree/master/batchRenaming_script/windows) and [Unix](https://github.com/jgarces02/FlowCT/tree/master/batchRenaming_script/linux)). It's comprised by two subscripts... see README(s) for additional information.

First step is unifying names for fluorofore channels and set them in the same order. This can be done through the *consolidate_and_reduction.rscript* script: it takes all FCS files inside a folder and make equal all channel names, puts them in the correct order and, if user says nothing, performs a events reduction. Running `Rscript consolidate_and_reduction.rscript --help` you can see default parameters and how to use them.

```
$ Rscript ../src/consolidate_and_reduction.rscript -h
Usage: consolidate_and_reduction.rscript [options]
Options:
        -d CHARACTER, --directory=CHARACTER
                path to FCS folder (default = current directory)
        -r LOGICAL, --reduction=LOGICAL
                logical indicating whether perform reduction (default = TRUE)
        -k INTEGER, --keep_n_events=INTEGER
                events to keep from original file (default = 10000)
        -s CHARACTER, --output_suffix=CHARACTER
                output suffix for new FCS files (default = '.ren')
        -o CHARACTER, --output_folder=CHARACTER
                output folder for storing new FCS files (default = current directory)
        -h, --help
                Show this help message and exit

$ Rscript ../src/consolidate_and_reduction.rscript
Processing file: MO_017294.fcs (keeping > 10000 events)
Processing file: MO_017564.fcs (keeping > 10000 events)
Processing file: MO_017612.fcs (keeping > 10000 events)
Processing file: MO_017809.fcs (keeping > 10000 events)
Processing file: MO_018207.fcs (keeping > 10000 events)
Processing file: MO_018208.fcs (keeping > 10000 events)
Processing file: SP_017294.fcs (keeping > 10000 events)
Processing file: SP_017565.fcs (keeping > 10000 events)
Processing file: SP_017612.fcs (keeping > 10000 events)
Processing file: SP_017809.fcs (keeping > 10000 events)
Processing file: SP_018207.fcs (keeping > 10000 events)
Processing file: SP_018208.fcs (keeping > 10000 events)

>>> Final markers order: FSC-A, FSC-H, SSC-A, SSC-H, FITC-A, PE-A, PerCP-Cy5-5-A, PE-Cy7-A, APC-A, APC-H7-A, V500-A, V450-A               
```

#### 0, bis) Setting R environment
```
setwd("~/FlowCT/data")

load_packages <- c("readxl", "flowCore", "flowAI", "flowViz", "flowStats", "gridExtra", "ggsci", "matrixStats", "ggplot2", "reshape2", 
                   "ggrepel", "dplyr", "RColorBrewer", "pheatmap", "FlowSOM", "ConsensusClusterPlus", "Rtsne", "uwot", "FlowCT",
                   "premessa", "phytools", "ggtree", "Hmisc", "corrplot", "ggthemes", "ggpubr", "matrixTests", "DataCombine")
invisible(lapply(load_packages, library, character.only = TRUE))
```

#### 1) Preprocessing and quality control
For this step is suposed that all your FCS have the same order and name for their fluorofore channels. If you haven't performed the reduction step previously, you can run: 
```
> reduce.FCS(file_or_directory = ".", keep_n_events = 10000)
FCS files within the selected folder:
[1] "SP_018207.fcs"     "SP_018208.fcs"     ...

Processing: SP_018207.fcs
Processing: SP_018208.fcs
...
```
In the case you've already run *consolidate_and_reduction.rscript* you can jump to next step.

Once we have events-reduced/names-checked/compensated files, we can start the quality check. To do that, we will use the package `flowAI`, an automatic method for the detection and removal of unwanted anomalies. Doublets removal is based on a ellipsoid gate in the FSC-A/FSC-H dotplot using the 95<sup>th</sup> percentile as limit for singlet cells (adjusted from Droumeva [one](https://github.com/LennonLab/flowcytometry/blob/master/bin/support_functions.R))
```
> qc.and.removeDoublets()
Processing: MO_017294.ren.red.fcs
Processing: MO_017564.ren.red.fcs
Processing: MO_017612.ren.red.fcs
Processing: MO_017809.ren.red.fcs
Processing: MO_018207.ren.red.fcs
Processing: MO_018208.ren.red.fcs
Processing: SP_017294.ren.red.fcs
Processing: SP_017565.ren.red.fcs
Processing: SP_017612.ren.red.fcs
Processing: SP_017809.ren.red.fcs
Processing: SP_018207.ren.red.fcs
Processing: SP_018208.ren.red.fcs
------------------------------
Final QC and remove doublets result:

|filenames | # initial events| % deleted events|Warning! |
|:---------|----------------:|----------------:|:--------|
|MO_017294 |            10000|             0.36|         |
|MO_017564 |            10000|             0.75|         |
|MO_017612 |            10000|             0.57|         |
|MO_017809 |            10000|             0.77|         |
|MO_018207 |            10000|             0.46|         |
|MO_018208 |            10000|             0.37|         |
|SP_017294 |            10000|             1.04|         |
|SP_017565 |            10000|             1.26|         |
|SP_017612 |            10000|             0.92|         |
|SP_017809 |            10000|             0.75|         |
|SP_018207 |            10000|             0.91|         |
|SP_018208 |            10000|             0.94|         |
```
During the process, if a sample losses more than 30% of events in the QC and the doublet removal, a *warning* (`(!)`) will appear in the final table. *Important, pay attention to suffixes from reduction step, if you use another different from default you must to adjust it in later steps this parameter (eg, `qc.and.removeDoublets(reduction_suffix = "_red"))`).*

#### 2) External compensation
#### 3) Importing, transformation and marker normalization
Here it's the reason for given an informative name to your FCS files: metadata generation is based on these filenames! We prefer a filename with a structure like *condition_patientID_somethingElse.fcs*. but of course, you can use your own code to generate this metadata or importing from an external table with additional information such as clinical data, for example. 
```
> (filenames <- list.files(pattern = ".preprocessed.fcs$"))
 [1] "MO_017294.ren.red.preprocessed.fcs" "MO_017564.ren.red.preprocessed.fcs"
 [3] "MO_017612.ren.red.preprocessed.fcs" "MO_017809.ren.red.preprocessed.fcs"
 [5] "MO_018207.ren.red.preprocessed.fcs" "MO_018208.ren.red.preprocessed.fcs"
 [7] "SP_017294.ren.red.preprocessed.fcs" "SP_017565.ren.red.preprocessed.fcs"
 [9] "SP_017612.ren.red.preprocessed.fcs" "SP_017809.ren.red.preprocessed.fcs"
[11] "SP_018207.ren.red.preprocessed.fcs" "SP_018208.ren.red.preprocessed.fcs"

> md <- data.frame(file_name = filenames, 
                 sample_id = 1:length(filenames),
                 condition = sapply(filenames,function(x) strsplit(x, split = "_|\\.")[[1]][1]),
                 patient_id = sapply(filenames,function(x) strsplit(x, split = "_|\\.")[[1]][2]))
> head(md)
                                                            file_name sample_id
MO_017294.ren.red.preprocessed.fcs MO_017294.ren.red.preprocessed.fcs         1
MO_017564.ren.red.preprocessed.fcs MO_017564.ren.red.preprocessed.fcs         2
MO_017612.ren.red.preprocessed.fcs MO_017612.ren.red.preprocessed.fcs         3
MO_017809.ren.red.preprocessed.fcs MO_017809.ren.red.preprocessed.fcs         4
MO_018207.ren.red.preprocessed.fcs MO_018207.ren.red.preprocessed.fcs         5
MO_018208.ren.red.preprocessed.fcs MO_018208.ren.red.preprocessed.fcs         6
                                   condition patient_id
MO_017294.ren.red.preprocessed.fcs        MO     017294
MO_017564.ren.red.preprocessed.fcs        MO     017564
MO_017612.ren.red.preprocessed.fcs        MO     017612
MO_017809.ren.red.preprocessed.fcs        MO     017809
MO_018207.ren.red.preprocessed.fcs        MO     018207
MO_018208.ren.red.preprocessed.fcs        MO     018208
```
And read all FCS into the same FlowSet objetc:

`fcs_raw <- read.flowSet(as.character(filenames), emptyValue = FALSE, transformation = FALSE, truncate_max_range = FALSE)`

Sometimes (more often we'd like :sweat:) markers names are different (eg, CCR4 and CD194). In with this function you can adjust them with the name you want:
```
>  markers.names(fcs_raw) #with no parameters, you only see them

|      |Florophore/Channel |Marker (default name) |New marker name  |
|:-----|:------------------|:---------------------|:----------------|
|$P1N  |FSC-A              |FSC-A                 |not defined yet! |
|$P2N  |FSC-H              |FSC-H                 |not defined yet! |
|$P3N  |SSC-A              |SSC-A                 |not defined yet! |
|$P4N  |SSC-H              |SSC-H                 |not defined yet! |
|$P5N  |FITC-A             |CD62L                 |not defined yet! |
|$P6N  |PE-A               |CXCR3                 |not defined yet! |
|$P7N  |PerCP-Cy5-5-A      |CD8                   |not defined yet! |
|$P8N  |CCR4:PE-Cy7-A      |CCR4                  |not defined yet! |
|$P9N  |APC-A              |CCR6                  |not defined yet! |
|$P10N |APC-H7-A           |CD4                   |not defined yet! |
|$P11N |V500-A             |CD45RA                |not defined yet! |
|$P12N |V450-A             |CD27                  |not defined yet! |

>  fcs_raw <- markers.names(fcs_raw, new_names = c("FSC_A", "FSC_H", "SSC_A", "SSC_H", "CD62L", "CXCR3", "CD8", "CD194", "CCR6 ", "CD4", "CD45", "CD27"))

|      |Florophore/Channel |Marker (default name) |New marker name |
|:-----|:------------------|:---------------------|:---------------|
|$P1N  |FSC-A              |FSC-A                 |FSC_A           |
|$P2N  |FSC-H              |FSC-H                 |FSC_H           |
|$P3N  |SSC-A              |SSC-A                 |SSC_A           |
|$P4N  |SSC-H              |SSC-H                 |SSC_H           |
|$P5N  |FITC-A             |CD62L                 |CD62L           |
|$P6N  |PE-A               |CXCR3                 |CXCR3           |
|$P7N  |PerCP-Cy5-5-A      |CD8                   |CD8             |
|$P8N  |CCR4:PE-Cy7-A      |CCR4                  |CD194           |
|$P9N  |APC-A              |CCR6                  |CCR6            |
|$P10N |APC-H7-A           |CD4                   |CD4             |
|$P11N |V500-A             |CD45RA                |CD45            |
|$P12N |V450-A             |CD27                  |CD27            |
```
And define those markers of interest and discard physical parameters for downstream analysis:
```
> surface_markers <- colnames(fcs_raw)[!grepl("FSC|SSC", colnames(fcs_raw))]
[1] "CD62L" "CXCR3" "CD8"   "CD194" "CCR6"  "CD4"   "CD45"  "CD27"
> H_markers <- colnames(fcs_raw)[grepl("FSC|SSC", colnames(fcs_raw))]
[1] "FSC_A" "FSC_H" "SSC_A" "SSC_H"
```
`FlowCore` allow us to transform flow cytometry data in a format that is more suitable for subsequent analysis. This transformation is necessary to deal with the bad representation of negative data values with digital cytometer: the Log transformation, indeed, leads to compression of data against the axes, poor visual representation of low intensity or unstained populations and errors in subsequent automatic gating/clustering algorithms... accordingly, we use the `arcsinh` transformation:
```
> cofact <- 500 #it should be optimized, if necessary, according to the experiment

> fcs <- fsApply(fcs_raw, function(x, cofactor = cofact){
  colnames(x)<- colnames(fcs_raw)
  expr <- exprs(x)
  expr <- asinh(expr[, c(surface_markers, H_markers)] / cofactor)
  exprs(x) <- expr
  x})

#generate a similar flowset without transformation for subsequent fcs generation
> fcs_no_transf<- fsApply(fcs_raw, function(x){
  colnames(x)<- colnames(fcs_raw)
  expr <- exprs(x)
  expr <- expr[, c(surface_markers, H_markers)] 
  exprs(x) <- expr
  x})
```
Despite the great efforts in the standardization of flow cytometry protocols, biological and technical variability can be reduced but not completely abrogated. To deal with this problem we used the `gaussNorm` method within the `flowStats` to reduce this technical variation by “shifting” and aligning the area of high local density (landmarks) in each channels. The first step in this process is to look at the distribution of peaks in each channel and identify parameters that needs to be aligned:
```
> densityplot(~ . , fcs, channels = c(surface_markers, H_markers), xlim = lims.FCS(fcs), filter = lapply(c(surface_markers, H_markers), curv1Filter))
```
![peaks_distro](https://github.com/jgarces02/FlowCT/blob/master/docs/normalization_first.png "Testing intensity peaks distro")

In our dataset it is easy to identify SSC_A and CD62L parameters as the ones with the highest variability between samples. After visual inspection you can decide to try normalization or exclude files from analysis, in our case we'll normalize CD62L and SSC_A:
```
> fcs_no_norm <- fcs #backup for non-normalized data
> markers_to_normalize <- c("SSC_A", "CD62L")

> for(marker in markers_to_normalize){
  datr <- gaussNorm(fcs, marker)$flowset
  if(require(flowViz)){
    grid.arrange(densityplot(as.formula(paste0("~", marker)), fcs, main = "original", xlim = lims.FCS(fcs), filter=curv1Filter(marker)),
                 densityplot(as.formula(paste0("~", marker)), datr, main = "normalized", xlim = lims.FCS(fcs), filter=curv1Filter(marker)),
                 ncol = 2)
  }
  fcs <- datr #storage normalized data
}
Adjusting the distance between landmarks
.....

Adjusting the distance between landmarks
.....
```
![peaks_distro_CD62L](https://github.com/jgarces02/FlowCT/blob/master/docs/CD62L_normalization.png "CD62L original and normalized peaks")

Finally, we can proceed with the generation of the “expression” matrix that would contain all the fluorescence intensity data at single cell level, and adjust initial metadata for each cell. 
```
> expr <- fsApply(fcs, exprs)
> expr_no_transf <- fsApply(fcs_no_transf, exprs)

#saturate marker intensities to values between 0 and 1... easier to interprete in heatmaps 
> rng <- colQuantiles(expr, probs = c(0.01, 0.99))
> expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
> expr01[expr01 < 0] <- 0
> expr01[expr01 > 1] <- 1

#generate sample IDs for each cell in the expression matrix and create a metadata table from initial one
> metadata_sc <- data.frame(sample_id = rep(md$sample_id, fsApply(fcs, nrow)), 
                          patient_id = rep(md$patient_id, fsApply(fcs, nrow)), 
                          condition = rep(md$condition, fsApply(fcs, nrow)))
> mdsc <- data.frame(metadata_sc, expr[,surface_markers])
> head(mdsc)
  sample_id patient_id condition     CD62L      CXCR3       CD8      CD194
1         1     017294        MO 1.3686484 0.13979632 0.1182595  1.6091200
2         1     017294        MO 0.1491844 1.45099442 4.1543690 -0.4079546
3         1     017294        MO 0.6164071 0.19323149 0.2415497  1.8109299
4         1     017294        MO 1.8904142 0.16306925 0.1473292  1.9848075
5         1     017294        MO 1.4291372 0.24543618 0.2426851  1.7537824
6         1     017294        MO 0.3724427 0.09374397 0.1159728  1.8338250
       CCR6         CD4      CD45       CD27
1 0.5761576 -0.09418564 0.3449939 0.44249329
2 0.9286080  0.41180576 0.2847564 3.55472430
3 0.8220239  0.97953895 0.6238420 0.07819405
4 0.8103202 -0.07624705 0.4794448 0.32365875
5 0.4633399 -0.37461350 0.2857165 0.43914630
6 0.7727667  0.57695024 0.3942684 0.29329506
```
If you have a very large amount of files it's almost impossible to visualize with the previously reported `densityplot` function, so a per marker overlay histograms could be very helpful:
```
> ggdf <- melt(mdsc, id.var = c("sample_id", "patient_id", "condition"), value.name = "expression", variable.name = "antigen")
> ggplot(ggdf, aes(x = expression, color = as.factor(sample_id))) +
    geom_density() + facet_wrap(. ~ antigen) + 
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 7), axis.text = element_text(size = 5))
```
![alt_viz_multiple_samples](https://github.com/jgarces02/FlowCT/blob/master/docs/global_hist.png "Alternative viz for multiple samples")

#### 4) Descriptive and exploratory analysis
Let's begin seeing the cell numbers for each condition:

`> cell.count.bx(data = metadata_sc, metadata = md)`

![boxplot_condition](https://github.com/jgarces02/FlowCT/blob/master/docs/cell_count_boxplot.metadata_sc.jpg "Boxplot by condition")

Dimensionality reduction algorithm aim to maintain the structure of the data (i.e. the distance between the different subpopulations, or patients in our case) while reducing the dimensions needed to observe differences. PCA is one of this techniques, and perform dimensionality reduction by projecting data to new coordinate preserving the variance. Mimicking a transcriptomic approach, Nowicka M et al8, proposed to use the median marker expression calculated over all cells. However, while this approach resulted very useful for Cytof data (and probably the same could be also for flow cytometry data with a high number of parameters), in our datasets the results were doubtful. Indeed, looking at the PCA as weell the heatmap in below figures it is easy to (wrongly) conclude that PB samples segregate differently from BM samples... we'll see later how using high number of events, ie. at single-cell level, results change.
```
#get the median marker expression per sample
> expr_median_sample_tbl <- data.frame(sample_id = metadata_sc$sample_id, expr) %>%
  group_by(sample_id) %>% summarize_all(list(median)) %>% as.data.frame()
> rownames(expr_median_sample_tbl) <- expr_median_sample_tbl$sample_id
> expr_median_sample_tbl <- expr_median_sample_tbl[,-1]

#compute the PCA
> drPCA <- dim.reduction(expr_data = expr_median_sample_tbl, metadata = md, reduction_method = "PCA")
> dr.plotting(drPCA$dr_melted, dr_calculated = "PCA", color_by = "condition", size = 3, output_type = "png", labels = "patient_id")

#draw heatmap
> annotation_col <- data.frame(row.names = rownames(expr_median_sample_tbl), condition = md$condition, patient_id = md$patient_id)
> color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100) #colors for the heatmap (expression values)
> x <- x2 <- c()
> set.seed(3); for(i in unique(md$condition)) x[i] <- sample(colors_palette, 1)
> set.seed(3); for(i in unique(md$patient_id)) x2[i] <- sample(colors_palette, 1)
> annotation_colors <- list(condition = x, patient_id = x2)

> pheatmap(t(expr_median_sample_tbl), color = color, display_numbers = FALSE,
         number_color = "black", fontsize_number = 5, annotation_col = annotation_col,
         annotation_colors = annotation_colors, clustering_method = "average")
```
<p float="left"> 
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/dr.PCA_col.condition.png" width="300" /> 
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/heatmap_patient.png" width="400" />  
</p> 

But changing at the single-cell level (previous subsampling to 1000 cells for each sample), PCA demonstrated an almost complete (and expected) overlapping between both cell distributions (PB and BM). Regarding the heatmap, we can observe a complete random distribution of cells derived from PB or BM, with the exception of a specific cluster of cells more abundant in the BM, presenting a complete negativity for the expression of all surface markers, thus probably being erythroblasts, obviously absent in PB.
# El subsampling se ha hecho con la función anterior, la actual da problemas... corregir
```
> sub_idx_heat <- sub.samples.idx(data = mdsc, colname_samples = "sample_id", samples_names = md$sample_id, subsampling = 1000, set.seed = 1234) #select how many cells to downsample per-sample
Extracting subsampling index for: 1
Extracting subsampling index for: 2
Extracting subsampling index for: 3
Extracting subsampling index for: 4
Extracting subsampling index for: 5
Extracting subsampling index for: 6
Extracting subsampling index for: 7
Extracting subsampling index for: 8
Extracting subsampling index for: 9
Extracting subsampling index for: 10
Extracting subsampling index for: 11
Extracting subsampling index for: 12
> heat_expr<- mdsc[sub_idx_heat,]

> ggdfPCA_1000 <- dim.reduction(expr_data = heat_expr[,4:ncol(heat_expr)], metadata = heat_expr[,1:3],  reduction_method = "PCA")
> dr.plotting(ggdfPCA_1000$dr_melted, dr_calculated = "PCA", color_by = "condition", output_type = NULL)
```
<p float="left"> 
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/scdr.PCA_col.condition.png" width="300" /> 
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/sc_heatmap_patient.png" width="400" />  
</p> 

#### 5) Cell clustering
As already stated before, and that's the reason for developin `FlowCT`, analyzing datasets of flow cytometry data with 8 or more colors (for exploratory/discovery purpose) with traditional gating strategies is absolutely not easy, time‐consuming and difficult to reproduce. Additionally, our capability to work with two-dimensional plots (or rarely with three-dimensional ones) only, renders quite impossible to keep an overview of all markers behavior on selected cell types, increasing the risk of missing valuable information. To overcome these limits, several methods of unsupervised cell clustering have been developed. A recent performance comparison of 18 clustering methods17 identified FlowSom18 as the first choice for exploratory analyses of large datasets, since this method was really fast and reproducible.
```
> fsom <- fsom.clustering(fcs, markers_to_use = "surface_markers", markers_to_plot = "tree", set.seed = 1234)
Calculating SOM clustering...
Building MST...
```
![MST_FlowSOM](https://github.com/jgarces02/FlowCT/blob/master/docs/MST_general.fcs.jpg "MST FlowSOM")

While `FlowSOM` already integrate the `metaClustering_consensus` function, to perform metaclustering we followed the approach suggested by Nowicka et al8 and used the ConsensusClusterPlus package to obtain a better control on each specific function. In this passage, we could decide the final number of clusters in which we would like our cells to be divided. It is important to select a number higher than the one expected: similar clusters could be later manually merged avoiding the risk to lose rare populations clusters. 

`> metaclusters <- fsom.metaclustering(fsom = fsom, num_clusters_metaclustering = 40, plotting = T, set.seed = 1234)`

![MST_metaclustering](https://github.com/jgarces02/FlowCT/blob/master/docs/MST_metaclustering.fsom.jpg "MST from metaclustering")

Below heatmap reporting the median values of each marker in each cluster can help (however in a limited way) to interprete results:

# Error en esta función, xq?: Error: Can't find columns `CD62L`, `CXCR3`, `CD8`, `CD194`, `CCR6`, … (and 3 more) in `.data`.
`> cluster.heatmap(expr = expr[, surface_markers], expr_saturated = expr01[, surface_markers],
                                cell_clusters = metaclusters$metaclusters)`

As previously described, dimensional reduction, togheter with SOM information, can be very useful to interprete cell behaviour giving us an idea of population distribution. Among the most recent algorithms for dimensionality reduction, t‑stochastic neighbour embedding (t-SNE)4,19 and Uniform Manifold Approximation and Projection (UMAP) represents the most popular ones. 
```
#combine single-cell data with SOM clusters
> mdsc_som <- data.frame(mdsc[,1:3], SOM = as.factor(metaclusters$metaclusters), mdsc[,4:ncol(mdsc)])

#subsampling
> sub_idx_som <- sub.samples.idx(data = mdsc_som, colname_samples = "sample_id", samples_names = md$sample_id, subsampling = 1000, set.seed = 1234)
> sel_expr <- mdsc_som[sub_idx_som,]

## Calculate dimensional reductions (simultaneously)
> dr <- dim.reduction(expr_data = sel_expr[,5:ncol(sel_expr)], metadata = sel_expr[,1:4], reduction_method = "all", set.seed = 1234)
Calculating PCA reduction...
Calculating tSNE reduction...
Calculating UMAP reduction...

#plot PCA, tSNE or UMAP ->  change "dr_calculated"
> dr.plotting(dr$dr_melted, dr_calculated = "tSNE", output_type = "png", color_by = "SOM")
> dr.plotting(dr$dr, dr_calculated = "UMAP", output_type = "png", color_by = "CD62L")

> for(facet in c("patient_id", "condition", "sample_id")){
  dr.plotting(dr$dr_melted, dr_calculated = "PCA", output_type = "png", color_by = "SOM", facet_by = facet, 
              output_name = paste0(facet, "_"))
}
PCA reduction colored by SOM and faceted by patient_id. Saved as -> patient_id_dr.PCA_col.SOM.png
PCA reduction colored by SOM and faceted by condition. Saved as -> condition_dr.PCA_col.SOM.png
PCA reduction colored by SOM and faceted by sample_id. Saved as -> sample_id_dr.PCA_col.SOM.png
```
<p float="left"> 
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/dr.tSNE_col.SOM.png" width="200" /> 
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/dr.UMAP_col.CD62L.png" width="200" />  
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/patient_id_dr.PCA_col.SOM.png" width="200" />  
</p> 

#### 6) Data exporting and external analysis
Once all the process of clustering and dimensional reduction have been completed, we can export our downsampled dataset in a new FCS files containing for each cell its expression value, dimensional reduction info (PCA, t-SNE and UMAP), SOM clusters... This file could be loaded in any flow cytometry software and subsequently analyzed to identify the population to which each cluster belongs.
```
> to_export <- data.frame(sample_id = as.numeric(sel_expr$sample_id), patient_id = as.numeric($patient_id), condition = as.numeric(sel_expr$condition), cluster_som = as.numeric(sel_expr$SOM),expr_no_transf[sub_idx_som,], dr$dr[,grepl("PCA|tSNE|UMAP", colnames(dr$dr))])

> outFile <- file.path(getwd(), "alltubeTh.fcs")
> write.FCS(as_flowFrame(as.matrix(to_export), source.frame = NULL), outFile)
[1] "~/FlowCT/data/results_preprocessing/alltubeTh.fcs"
```
..... external tutorial? ......

The final step for external analysis is to create an excel file (or a text file and modify consequently code below) assigning a new name for each cluster, some like that...

|original_cluster | new_cluster     |
|----------------:|:----------------|
|1                | debris          |
|2                | lymphocytes     |
|3                | lymphocytes     |
|4                | eosinophils     |
|5                | erythroblasts   |
|...              | ...             |

...read it in R and replace previuous clusters with new ones, drawing again dimensional reduction plots:
```
> cluster_merging1 <- read_excel("../Th_new_clustering.xlsx", col_types = "text")

> cell_clustering1m <- FindReplace(data.frame(cell_clustering1m = as.factor(metaclusters$metaclusters)), Var = "cell_clustering1m", replaceData = cluster_merging1, from = "original_cluster", to = "new_cluster", vector = T) %>% as.factor()
Only exact matches will be replaced.

> cell_clustering1m_plotStars <- FindReplace(data.frame(cell_clustering1m = as.factor(metaclusters$plotStars_value)), Var = "cell_clustering1m", replaceData = cluster_merging1, from = "original_cluster", to = "new_cluster", vector = T) %>% as.factor()
Only exact matches will be replaced.

> dr$dr_melted$cell_clustering1m <- FindReplace(data.frame(cell_clustering1m = as.factor(dr$dr$SOM)), 
                                    Var = "cell_clustering1m", replaceData = cluster_merging1, from = "original_cluster", to = "new_cluster", vector = T) %>% as.factor()
Only exact matches will be replaced.

> dr.plotting(dr$dr_melted, dr_calculated = "tSNE", color_by = "cell_clustering1m", output_type = "png")
tSNE reduction colored by cell_clustering1m. Saved as -> dr.tSNE_col.cell_clustering1m.png
> dr.plotting(dr$dr_melted, dr_calculated = "UMAP", color_by = "cell_clustering1m", output_type = "png")
UMAP reduction colored by cell_clustering1m. Saved as -> dr.UMAP_col.cell_clustering1m.png

> PlotStars(fsom, backgroundValues = cell_clustering1m_plotStars, backgroundColor = alpha(colors_palette, alpha = 0.4))

> cluster_heatmap(expr = expr[, surface_markers], expr_saturated = expr01[, surface_markers], cell_clusters = cell_clustering1m)
```
<p float="left"> 
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/dr.tSNE_col.cell_clustering1m.png" width="200" /> 
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/dr.UMAP_col.cell_clustering1m.png" width="200" />  
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/MST_newClust.png" width="200" />  
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/heatmap......" width="200" />  
</p> 

#### 7) Subclustering analysis 
At this point, the next step is to perform a deeper analysis on specific subpopulation identified with the FlowSom clustering followed by cluster merging... so, let's go to repeat the above analysis (very quickly)!.
```
#generation of a new matrix for subclustering
> metadata_sc <- cbind(mdsc, FlowSOM = cell_clustering1m)

#keep only desired cell population
> matrix_clusters <- as.matrix(metadata_sc[metadata_sc$FlowSOM == "lymphocytes", surface_markers])

> rngL <- colQuantiles(matrix_clusters, probs = c(0.01, 0.99))
> matrix_clusters01 <- t((t(matrix_clusters) - rngL[, 1]) / (rngL[, 2] - rngL[, 1]))
> matrix_clusters01[matrix_clusters01 < 0] <- 0
> matrix_clusters01[matrix_clusters01 > 1] <- 1

> metadata_scL <- metadata_sc[metadata_sc$FlowSOM == "lymphocytes",]
> rownames(metadata_scL) <- NULL #reset rownames after cell selection (for later subsampling)

> expr_no_transfL <- expr_no_transf[metadata_sc$FlowSOM == "lymphocytes",]

#lymphocytes subclustering
> cell.count.bx(metadata_scL, counts_by = "sample_id", metadata = md)
```
![boxplot_subclust](https://github.com/jgarces02/FlowCT/blob/master/docs/MST_metaclustering.fsom.jpg "Boxplot specific subclustering")
```
#FlowSOM
> fsomL <- fsom.clustering(data = metadata_scL[,surface_markers], markers_to_use = "surface_markers", set.seed = 1234)
> metaclustersL <- fsom.metaclustering(fsom = fsomL, num_clusters_metaclustering = 40, plotting = T, set.seed = 1234)

> cluster_heatmap(expr = metadata_scL[,surface_markers], expr_saturated = matrix_clusters01[,surface_markers], cell_clusters = metaclustersL$metaclusters)

> metadata_scL$FlowSOM_L <- metaclustersL$metaclusters
```
imagen...
```
#dim reduction
> sub_idxL <- sub.samples.idx(metadata_scL, colname_samples = "sample_id", samples_names = md$sample_id, subsampling = 1000, set.seed = 1234)
> sel_exprL <- metadata_scL[sub_idxL,]

> drL <- dim.reduction(sel_exprL[,surface_markers], metadata = sel_exprL[,!(colnames(sel_exprL) %in% surface_markers)], reduction_method = "all", set.seed = 1234)

> gg_list <- list() #list to store plots
> for(marker in surface_markers) gg_list[[marker]] <- dr.plotting(drL$dr, dr_calculated = "PCA", color_by = marker, output_type = NULL)
> do.call("grid.arrange", c(gg_list, nrow = 2)) #show all plots
```
![dr_markers_subclust](https://github.com/jgarces02/FlowCT/blob/master/docs/marker_expression_sublust.png "Multimarker ploting dr subclustering")
# Algo está petando en la sustitución en el DR...
```
#exporting
> to_export <- data.frame(sample_id = as.numeric(sel_exprL$sample_id), patient_id = as.numeric(sel_exprL$patient_id), condition = as.numeric(sel_exprL$condition), cluster_somL = as.numeric(sel_exprL$FlowSOM_L), expr_no_transfL[sub_idxL,], drL$dr[,grepl("PCA|tSNE|UMAP", colnames(drL$dr))])
> outFileL <- file.path(getwd(), "alltubeThclustL.fcs")
> write.FCS(as_flowFrame(as.matrix(to_export), source.frame = NULL), outFileL)

# ++++++++++++++++++++++++++++
# external step >>> manual analysis in a flow cytometry software to identify clusters
# ++++++++++++++++++++++++++++

#merge clusters
> cluster_merging1L <- read_excel("../../../Th_new_clusteringL.xlsx", col_types = "text")

> cell_clustering1mL <- FindReplace(data.frame(cell_clustering1m = as.factor(metaclustersL$metaclusters)), Var = "cell_clustering1m", replaceData = cluster_merging1L, from = "original_cluster", to = "new_cluster", vector = T) %>% as.factor()
> cell_clustering1mL_plotStars <- FindReplace(data.frame(cell_clustering1m = as.factor(metaclustersL$plotStars_value)), Var = "cell_clustering1m", replaceData = cluster_merging1L, from = "original_cluster", to = "new_cluster", vector = T) %>% as.factor()
> drL$dr_melted$cell_clustering1mL <- FindReplace(data.frame(cell_clustering1m = as.factor(drL$dr$FlowSOM_L)), Var = "cell_clustering1m", replaceData = cluster_merging1L, from = "original_cluster", to = "new_cluster", vector = T) %>% as.factor()

#new dr plots
> PlotStars(fsomL, backgroundValues = cell_clustering1mL_plotStars,  backgroundColor = alpha(colors_palette, alpha = 0.5))

> dr.plotting(drL$dr_melted, dr_calculated = "tSNE", color_by = "cell_clustering1mL", output_type = NULL, facet_by = "patient_id")

> cluster_heatmap(expr = matrix_clusters[,surface_markers], expr_saturated = matrix_clusters01[,surface_markers], cell_clusters = cell_clustering1mL)
```
To have a better idea of distribution and phenotype of each population we represented our results in a phylogenic tree taking advantage of the recently developed ggtree22 package. The heatmap associated to each cell population reports the median expression of each marker, while the dimension of the circle is proportional to the relative abundance of the population itself. Next we performed a single cell heatmap to verify if the manual cluster merging respected the unsupervised distribution of cells.

First line (`circ.tree.selectNodes`) draws a simple circular dendrogram to select those nodes where we can cut and differentiate each cell subpopulation (that will be colored later with `circ.tree`).
```
> circ.tree.selectNodes(exprs = matrix_clusters, exprs_saturated = matrix_clusters01, cell_clusters = cell_clustering1mL)
> circ.tree(exprs = matrix_clusters, exprs_saturated = matrix_clusters01, cell_clusters = cell_clustering1mL, dendro_labels = F, nodes = c(29,33,35))

#heatmap
> sub_idxL_heat <- sub.samples.idx(metadata_scL, colname_samples = "sample_id", subsampling = 1000, set.seed = 1234)
> sel_exprL_heat <- metadata_scL[sub_idxL_heat,]

> annotation_col <- data.frame(row.names = rownames(sel_exprL_heat), sel_exprL_heat[,!(colnames(sel_exprL_heat) %in% surface_markers)])
> annotation_colors <- list(condition = c(BM = colors_palette[1], PB = colors_palette[2]))

> pheatmap(t(sel_exprL_heat[,surface_markers]), annotation_col = annotation_col, clustering_method = "average", color = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100), display_numbers = FALSE, number_color = "black", fontsize_number = 5, annotation_colors = annotation_colors, show_colnames = F, treeheight_col = 0, treeheight_row = 0)
```
