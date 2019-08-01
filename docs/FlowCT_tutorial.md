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
4. Descriptive analysis through heatmaps, boxplots, PCA... and general and single-cell level.
5. Cell clustering using SOM (Self-Organizing Map), PCA, tSNE and UMAP and respective visualizations.
6. New global FCS generation and manually identification of generated clusters in the previous step. Later, assign and merge of identified cell populations.
7. (Optional) Subclustering and repeeat the analysis on a desired specific population.
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
>reduce.FCS(file_or_directory = ".", keep_n_events = 10000)
FCS files within the selected folder:
[1] "SP_018207.fcs"     "SP_018208.fcs"     ...

Processing: SP_018207.fcs
Processing: SP_018208.fcs
...
```
In the case you've already run *consolidate_and_reduction.rscript* you can jump to next step.

Once we have events-reduced/names-checked/compensated files, we can start the quality check. To do that, we will use the package `flowAI`, an automatic method for the detection and removal of unwanted anomalies. Doublets removal is based on a ellipsoid gate in the FSC-A/FSC-H dotplot using the 95<sup>th</sup> percentile as limit for singlet cells (adjusted from Droumeva [one](https://github.com/LennonLab/flowcytometry/blob/master/bin/support_functions.R))
```
>qc.and.removeDoublets()
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
>(filenames <- list.files(pattern = ".preprocessed.fcs$"))
 [1] "MO_017294.ren.red.preprocessed.fcs" "MO_017564.ren.red.preprocessed.fcs"
 [3] "MO_017612.ren.red.preprocessed.fcs" "MO_017809.ren.red.preprocessed.fcs"
 [5] "MO_018207.ren.red.preprocessed.fcs" "MO_018208.ren.red.preprocessed.fcs"
 [7] "SP_017294.ren.red.preprocessed.fcs" "SP_017565.ren.red.preprocessed.fcs"
 [9] "SP_017612.ren.red.preprocessed.fcs" "SP_017809.ren.red.preprocessed.fcs"
[11] "SP_018207.ren.red.preprocessed.fcs" "SP_018208.ren.red.preprocessed.fcs"

>md <- data.frame(file_name = filenames, 
                 sample_id = 1:length(filenames),
                 condition = sapply(filenames,function(x) strsplit(x, split = "_|\\.")[[1]][1]),
                 patient_id = sapply(filenames,function(x) strsplit(x, split = "_|\\.")[[1]][2]))
>head(md)
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
> markers.names(fcs_raw) #with no parameters, you only see them

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

> fcs_raw <- markers.names(fcs_raw, new_names = c("FSC_A", "FSC_H", "SSC_A", "SSC_H", "CD62L", "CXCR3", "CD8", "CD194", "CCR6 ", "CD4", "CD45", "CD27"))

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
>surface_markers <- colnames(fcs_raw)[!grepl("FSC|SSC", colnames(fcs_raw))]
[1] "CD62L" "CXCR3" "CD8"   "CD194" "CCR6"  "CD4"   "CD45"  "CD27"
>H_markers <- colnames(fcs_raw)[grepl("FSC|SSC", colnames(fcs_raw))]
[1] "FSC_A" "FSC_H" "SSC_A" "SSC_H"
```
`FlowCore` allow us to transform flow cytometry data in a format that is more suitable for subsequent analysis. This transformation is necessary to deal with the bad representation of negative data values with digital cytometer: the Log transformation, indeed, leads to compression of data against the axes, poor visual representation of low intensity or unstained populations and errors in subsequent automatic gating/clustering algorithms... accordingly, we use the `arcsinh` transformation:
```
>cofact <- 500 #it should be optimized, if necessary, according to the experiment

>fcs <- fsApply(fcs_raw, function(x, cofactor = cofact){
  colnames(x)<- colnames(fcs_raw)
  expr <- exprs(x)
  expr <- asinh(expr[, c(surface_markers, H_markers)] / cofactor)
  exprs(x) <- expr
  x})

#generate a similar flowset without transformation for subsequent fcs generation
>fcs_no_transf<- fsApply(fcs_raw, function(x){
  colnames(x)<- colnames(fcs_raw)
  expr <- exprs(x)
  expr <- expr[, c(surface_markers, H_markers)] 
  exprs(x) <- expr
  x})
```
Despite the great efforts in the standardization of flow cytometry protocols, biological and technical variability can be reduced but not completely abrogated. To deal with this problem we used the `gaussNorm` method within the `flowStats` to reduce this technical variation by “shifting” and aligning the area of high local density (landmarks) in each channels. The first step in this process is to look at the distribution of peaks in each channel and identify parameters that needs to be aligned:
```
>densityplot(~ . , fcs, channels = c(surface_markers, H_markers), xlim = lims.FCS(fcs), filter = lapply(c(surface_markers, H_markers), curv1Filter))
```
![peaks_distro](https://github.com/jgarces02/FlowCT/blob/master/docs/normalization_first.png "Testing intensity peaks distro")

In our dataset it is easy to identify SSC_A and CD62L parameters as the ones with the highest variability between samples. After visual inspection you can decide to try normalization or exclude files from analysis, in our case we'll normalize CD62L and SSC_A:
```
>fcs_no_norm <- fcs #backup for non-normalized data
>markers_to_normalize <- c("SSC_A", "CD62L")

>for(marker in markers_to_normalize){
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
>expr <- fsApply(fcs, exprs)
>expr_no_transf <- fsApply(fcs_no_transf, exprs)

#saturate marker intensities to values between 0 and 1... easier to interprete in heatmaps 
>rng <- colQuantiles(expr, probs = c(0.01, 0.99))
>expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
>expr01[expr01 < 0] <- 0
>expr01[expr01 > 1] <- 1

#generate sample IDs for each cell in the expression matrix and create a metadata table from initial one
>metadata_sc <- data.frame(sample_id = rep(md$sample_id, fsApply(fcs, nrow)), 
                          patient_id = rep(md$patient_id, fsApply(fcs, nrow)), 
                          condition = rep(md$condition, fsApply(fcs, nrow)))
>mdsc <- data.frame(metadata_sc, expr[,surface_markers])
>head(mdsc)
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
>ggdf <- melt(mdsc, id.var = c("sample_id", "patient_id", "condition"), value.name = "expression", variable.name = "antigen")
>ggplot(ggdf, aes(x = expression, color = as.factor(sample_id))) +
    geom_density() + facet_wrap(. ~ antigen) + 
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 7), axis.text = element_text(size = 5))
```
![alt_viz_multiple_samples](https://github.com/jgarces02/FlowCT/blob/master/docs/global_hist.png "Alternative viz for multiple samples")

#### 4) Descriptive and exploratory analysis
Let's begin seeing the cell numbers for each condition:

`cell.count.bx(data = metadata_sc, metadata = md)`

![boxplot_condition](https://github.com/jgarces02/FlowCT/blob/master/docs/cell_count_boxplot.metadata_sc.jpg "Boxplot by condition")

Dimensionality reduction algorithm aim to maintain the structure of the data (i.e. the distance between the different subpopulations, or patients in our case) while reducing the dimensions needed to observe differences. PCA is one of this techniques, and perform dimensionality reduction by projecting data to new coordinate preserving the variance. Mimicking a transcriptomic approach, Nowicka M et al8, proposed to use the median marker expression calculated over all cells. However, while this approach resulted very useful for Cytof data (and probably the same could be also for flow cytometry data with a high number of parameters), in our datasets the results were doubtful. Indeed, looking at the PCA as weell the heatmap in below figures it is easy to (wrongly) conclude that PB samples segregate differently from BM samples... we'll see later how using high number of events, ie. at single-cell level, results change.
```
#get the median marker expression per sample
>expr_median_sample_tbl <- data.frame(sample_id = metadata_sc$sample_id, expr) %>%
  group_by(sample_id) %>% summarize_all(list(median)) %>% as.data.frame()
>rownames(expr_median_sample_tbl) <- expr_median_sample_tbl$sample_id
>expr_median_sample_tbl <- expr_median_sample_tbl[,-1]

#compute the PCA
>drPCA <- dim.reduction(expr_data = expr_median_sample_tbl, metadata = md, reduction_method = "PCA")
>dr.plotting(drPCA$dr_melted, dr_calculated = "PCA", color_by = "condition", size = 3, output_type = "png", labels = "patient_id")

#draw heatmap
>annotation_col <- data.frame(row.names = rownames(expr_median_sample_tbl), condition = md$condition, patient_id = md$patient_id)
>color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100) #colors for the heatmap (expression values)
>x <- x2 <- c()
>set.seed(3); for(i in unique(md$condition)) x[i] <- sample(colors_palette, 1)
>set.seed(3); for(i in unique(md$patient_id)) x2[i] <- sample(colors_palette, 1)
>annotation_colors <- list(condition = x, patient_id = x2)

>pheatmap(t(expr_median_sample_tbl), color = color, display_numbers = FALSE,
         number_color = "black", fontsize_number = 5, annotation_col = annotation_col,
         annotation_colors = annotation_colors, clustering_method = "average")
```
<p float="left">
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/scdr.PCA_col.condition.png" width="400" />
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/sc_heatmap_patient.png" width="400" /> 
</p>

But changing at the single-cell level (previous subsampling to 1000 cells for each sample), PCA demonstrated an almost complete (and expected) overlapping between both cell distributions (PB and BM). Regarding the heatmap, we can observe a complete random distribution of cells derived from PB or BM, with the exception of a specific cluster of cells more abundant in the BM, presenting a complete negativity for the expression of all surface markers, thus probably being erythroblasts, obviously absent in PB.
# El subsampling se ha hecho con la función anterior, la actual da problemas... corregir
```
>sub_idx_heat <- sub.samples.idx(data = mdsc, colname_samples = "sample_id", samples_names = md$sample_id, subsampling = 1000, set.seed = 1234) #select how many cells to downsample per-sample
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
>heat_expr<- mdsc[sub_idx_heat,]

>ggdfPCA_1000 <- dim.reduction(expr_data = heat_expr[,4:ncol(heat_expr)], metadata = heat_expr[,1:3],  reduction_method = "PCA")
>dr.plotting(ggdfPCA_1000$dr_melted, dr_calculated = "PCA", color_by = "condition", output_type = NULL)
```
<p float="left">
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/scdr.PCA_col.condition.png" width="400" />
  <img src="https://github.com/jgarces02/FlowCT/blob/master/docs/sc_heatmap_patient.png" width="400" /> 
</p>
