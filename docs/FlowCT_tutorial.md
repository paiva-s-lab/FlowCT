## FlowCT (R package)
### A semi-automated workflow for deconvolution of immunophenotypic data and objective reporting on large datasets

Hello! Thanks for interesting you in FlowCT!

This is a very generic tutorial showing some possibilites of this R package. For any additional information or check another analysis appoximations, you can see our paper in [xxxxx](wwww.web.com) and the help for each function of this package.

# Adjust this part with general pipeline below!
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
```

#### 0, bis) Setting R environment
```
setwd("~/projects/FlowCT_tutorial")

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
[1] "SP_018207.fcs"     "SP_018208.fcs"

Processing: SP_018207.fcs
Processing: SP_018208.fcs
```
In the case you've already run *consolidate_and_reduction.rscript* you can jump to next step.

Once we have events-reduced/names-checked/compensated files, we can start the quality check. To do that, we will use the package `flowAI`, an automatic method for the detection and removal of unwanted anomalies. Doublets removal is based on a ellipsoid gate in the FSC-A/FSC-H dotplot using the 95<sup>th</sup> percentile as limit for singlet cells (adjusted from Droumeva [one](https://github.com/LennonLab/flowcytometry/blob/master/bin/support_functions.R))
```
>qc.and.removeDoublets()
Processing: SP_018207.red.fcs
Processing: SP_018208.red.fcs
------------------------------
Final QC and remove doublets result:

|filenames     | # initial events| % deleted events|Warning! |
|:-------------|----------------:|----------------:|:--------|
|SP_018207.red |            10000|             0.96|         |
|SP_018208.red |            10000|             1.12|         |
```
During the process, if a sample losses more than 30% of events in the QC and the doublet removal, a *warning* (`(!)`) will appear in the final table. *Important, pay attention to suffixes from reduction step, if you use another different from default you must to adjust it in later steps this parameter (eg, `qc.and.removeDoublets(reduction_suffix = "_red"))`).*

#### 2) External compensation
#### 3) Importing, transformation and marker normalization
Here it's the reason for given an informative name to your FCS files: metadata generation is based on these filenames! We prefer a filename with a structure like *condition_patientID_somethingElse.fcs*. but of course, you can use your own code to generate this metadata or importing from an external table with additional information such as clinical data, for example. 
```
>(filenames <- list.files(pattern = ".preprocessed.fcs$"))
[1] "SP_018207_red.preprocessed.fcs" "SP_018208_red.preprocessed.fcs"

>md <- data.frame(file_name = filenames, 
                 sample_id = 1:length(filenames),
                 condition = sapply(filenames,function(x) strsplit(x, split = "_|\\.")[[1]][1]),
                 patient_id = sapply(filenames,function(x) strsplit(x, split = "_|\\.")[[1]][2]))
>head(md)
                                                    file_name sample_id
SP_018207_red.preprocessed.fcs SP_018207_red.preprocessed.fcs         1
SP_018208_red.preprocessed.fcs SP_018208_red.preprocessed.fcs         2
                               condition patient_id
SP_018207_red.preprocessed.fcs        SP     018207
SP_018208_red.preprocessed.fcs        SP     018208
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

In our dataset it is easy to identify ¿SSC_A? and CD62L parameters as the ones with the highest variability between samples. After visual inspection you can decide to try normalization or exclude files from analysis, in our case we'll normalize CD62L and ¿SSC_A?:
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
  sample_id patient_id condition        CD62L      CXCR3         CD8      CD194
1         1     018207        SP  0.100415325 0.12565532  0.09410983 2.88623478
2         1     018207        SP  1.684631632 0.21251500  0.33678475 2.43242198
3         1     018207        SP  2.713427007 0.21018863 -0.08113977 1.41647705
4         1     018207        SP  0.106943000 0.08280282 -0.12873575 0.07693925
5         1     018207        SP -0.001977132 0.35285651 -0.13759982 3.30589471
6         1     018207        SP  0.525919523 0.06810811  0.37944536 0.43635330
       CCR6         CD4       CD45      CD27
1 0.1410125 -0.10364878 0.02049048 0.3859293
2 0.6602614  2.75647495 1.01969186 0.3241081
3 0.9508190  0.20926000 0.59292335 0.2593487
4 0.0492510 -0.01760147 0.51407507 1.1442070
5 0.5231158  0.08819571 0.13633504 0.2543581
6 0.6266562  0.30556337 0.56742354 0.3538544
```
