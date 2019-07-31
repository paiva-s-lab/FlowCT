## FlowCT (R package)
### A semi-automated workflow for deconvolution of immunophenotypic data and objective reporting on large datasets

Hello! Thanks for interesting you in FlowCT!

This is a very generic tutorial showing some possibilites of this R package. For any additional information or check another analysis appoximations, you can see our paper in [xxxxx](wwww.web.com) and the help for each function of this package.

So let's go! Pipeline is structurated in multiple parts, some internals and other externals:
1. Preprocess of FCS files, ie. to unify (florofore) channels names and set them in the same order. Normaly, it depends of computation power available, user also wants to reduce the number of events to avoid computer crashing or to facilitate compensation of multiple FCS.
2. External compensation of FCS files (if needed) to correct overlap between channels. Whatever flow cytometry software can be used, in our case we use [Infinicyt<sup>TM</sup>](http://www.infinicyt.com/).
3. Import data from all FCS compensated files, check uniformity of marker names (eg. CCR4 and not CD194), quality control, doublet removal and transformation/normalization.
4. Descriptive analysis through heatmaps, boxplots, PCA... and general and single-cell level.
5. Cell clustering using SOM (Self-Organizing Map), PCA, tSNE and UMAP and respective visualizations.
6. New global FCS generation and manually identification of generated clusters in the previous step. Later, assign and merge of identified cell populations.
7. (Optional) Subclustering and repeeat the analysis on a desired specific population.
8. Exporting table for statistical analysis, either inside R or another statistic software.

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

#### 1) Preprocessing
For this step is suposed that all your FCS have the same order and name for their fluorofore channels. If you haven't performed the reduction step previously, you can run: 
```
>reduce.FCS(file_or_directory = ".", keep_n_events = 10000)
FCS files within the selected folder:
[1] "SP_018207.fcs"     "SP_018208.fcs"

Processing: SP_018207.fcs
Processing: SP_018208.fcs
```
Following it's the quality control and the removal of doublets. If a sample losses more than 30% of events in the QC and the doublet removal, a *warning* (`(!)`) will appear in the final table.
```
>qc.and.removeDoublets()
Processing: SP_018207_red.fcs
Processing: SP_018208_red.fcs
------------------------------
Final QC and remove doublets result:

|filenames     | # initial events| % deleted events|Warning! |
|:-------------|----------------:|----------------:|:--------|
|SP_018207_red |            10000|             0.96|         |
|SP_018208_red |            10000|             1.12|         |
```

