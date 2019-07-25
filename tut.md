### Backbone for tutorial writting...


```
#################################################
### Environment setting #########################
#################################################
## Working directory
setwd("G:/Mi unidad/Proyectos/FlowCT_Ciro/")

## Load packages and functions
load_packages <- c("readxl", "flowCore", "flowAI", "flowViz", "flowStats", "gridExtra", "ggsci", "matrixStats", "ggplot2", "reshape2", 
                   "ggrepel", "dplyr", "RColorBrewer", "pheatmap", "FlowSOM", "ConsensusClusterPlus", "Rtsne", "uwot", 
                   "premessa", "phytools", "ggtree", "Hmisc", "corrplot", "ggthemes", "ggpubr", "matrixTests", "DataCombine")
lapply(load_packages, suppressPackageStartupMessages(library), character.only = TRUE)

# source("functions_FlowCT_v2.4.R")
# options(rsconnect.http = "internal")
# Sys.setenv(http_proxy  = "http://proxy.unav.es:8080")
# Sys.setenv(https_proxy = "http://proxy.unav.es:8080")
devtools::install_github("jgarces02/p", auth_token = "21ea9880f944d42755479e54a5b19ddd00fe17f6")
library(FlowCT)

################################################
### Pre-processing of FCS data files ############
#################################################

## Standarize names for markers (GUI)
# paneleditor_GUI() #it's not working properly yet...

## Downsample FCS files (for computer requirements) and standarize CD names
# setwd("files_FlowCT/renamed/")
setwd("files_FlowCT_new/")
reduce.FCS(file_or_directory = ".", keep_n_events = 10000)

## Quality control and doublet removal
setwd("reduced/")
qc.and.removeDoublets(reduction_computed = T)

setwd("results_HQsinglets/")

##################################################
### Generation of flowset objects to be analyzed #
##################################################
## Metadata matrix construction
# >> NOTE: this is performed using the filename of each file, it must be named following this structure: condition_patientID_somethingElse.fcs
(filenames <- list.files(pattern = "_preprocessed.fcs$"))

md <- data.frame(file_name = filenames, 
                 sample_id = 1:length(filenames),
                 condition = sapply(filenames,function(x) strsplit(x, split = "_")[[1]][1]),
                 patient_id = sapply(filenames,function(x) strsplit(x, split = "_")[[1]][2]))
head(md)

## Flowset generation
fcs_raw <- read.flowSet(as.character(md$file_name), emptyValue = FALSE, transformation = FALSE, truncate_max_range = FALSE)

#adjust maker names
markers.names(fcs_raw)
fcs_raw <- markers.names(fcs_raw, 
                         new_names = c("FSC_A", "FSC_H", "SSC_A", "SSC_H", "CD62L", "CXCR3", "CD8", "CCR4", "CCR6", "CD4", "CD27", "CD45RA"))


#read markers panel used
# panel <- read_excel("/mnt/beegfs/cdalisvam/paper/Th_panel.xlsx") 
panel <- read_excel("../../Th_panel.xlsx") 
head(data.frame(panel))

panel$Antigen <- gsub("-", "_", panel$Antigen) #replace problematic characters

surface_markers <- panel$Antigen[panel$use == 1] #select markers of interest
H_markers <- panel$Antigen[panel$use == 0]

#check if panel markers imported from Excel (ie, panel) are equal to read FCSs (ie, fcs_raw)
markers.equal(surface_markers, fcs_raw@colnames)

## Flowset transformation (arcsinh)
cofact <- 500 #default. It should be optimized, if necessary, according to the experiment

fcs <- fsApply(fcs_raw, function(x, cofactor = cofact){
  colnames(x)<- colnames(fcs_raw)
  expr <- exprs(x)
  expr <- asinh(expr[, c(surface_markers, H_markers)] / cofactor)
  exprs(x) <- expr
  x})

#generate a similar flowset without transformation for subsequent fcs generation
fcs_no_transf<- fsApply(fcs_raw, function(x){
  colnames(x)<- colnames(fcs_raw)
  expr <- exprs(x)
  expr <- expr[, c(surface_markers, H_markers)] 
  exprs(x) <- expr
  x})

# save(list = c("fcs", "fcs_no_transf"), file = "raw_and_transformed_data.Rdata") #save raw and transformed data

## Normalization and data alignment
# is normalization needed? first evaluation of histograms to decide on it!
# for high number of samples use represetation of line 202 (overlay histograms)

q <- c("CD62L", "CXCR3","CD8","CCR4","CCR6","CD4","CD27","CD45RA")
densityplot(~ ., fcs, channels = q, xlim = lims.FCS(fcs),
            filter = lapply(q, curv1Filter))


#after visual inspection you can decide to try to normalize or exclude files from analysis!
fcs_no_norm <- fcs #backup for non-normalized data
markers_to_normalize <- c("SSC_A", "CD62L")

for(marker in markers_to_normalize){
  datr <- gaussNorm(fcs, marker)$flowset
  if(require(flowViz)){
    grid.arrange(densityplot(as.formula(paste0("~", marker)), fcs, main = "original", xlim = lims.FCS(fcs), filter=curv1Filter(marker)),
                 densityplot(as.formula(paste0("~", marker)), datr, main = "normalized", xlim = lims.FCS(fcs), filter=curv1Filter(marker)),
                 ncol = 2)
  }
  fcs <- datr #storage normalized data
}

## Extraction of single-cell fluorescence values and generation of expression matrixes
#extract single cell values for expression
expr <- fsApply(fcs, exprs)
head(expr)

expr_no_transf <- fsApply(fcs_no_transf, exprs)
head(expr_no_transf)

## According to......we decided to follow the suggestion to apply a second transformation that scales expression of all 
##markers to values between 0 and 1. This make values easier to be analyzed in heatmaps 
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

#generate sample IDs for each cell in the expression matrix and create a metadata table
metadata_sc <- data.frame(sample_id = rep(md$sample_id, fsApply(fcs, nrow)), 
                          patient_id = rep(md$patient_id, fsApply(fcs, nrow)), 
                          condition = rep(md$condition, fsApply(fcs, nrow)))
head(metadata_sc)

#combine metadata and expression data at single cell level
mdsc <- data.frame(metadata_sc, expr[,surface_markers])
head(mdsc)

#################################################
### Descriptive and exploratory analysis ########
#################################################
## Overlay histograms and total cell count per group
ggdf <- melt(mdsc, id.var = c("sample_id", "patient_id", "condition"), value.name = "expression", variable.name = "antigen")
head(ggdf)

markers_normalized <- c("SSC_A", "CD62L")
for(marker in markers_normalized){
  print(ggplot(ggdf[ggdf$antigen == marker,], aes(x = expression, color = condition, group = sample_id)) + #change "data" if makers individually
    geom_density() +
    theme_minimal() + ggtitle(marker) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
    scale_color_manual(values = colors_palette))
}

# ggplot(data = ggdf, aes(x = expression, color = condition, group = sample_id)) +  #plot all markers at the same time
#   geom_density(size = 1.5) +
#   facet_wrap(~ antigen, nrow = 4, scales = "free") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1),
#         strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
#   scale_color_manual(values = colors_palette)

#boxplot showing the number of cells for condition
cell.count.bx(data = metadata_sc, metadata = md)

## PCA and heatmap graphical representations
#get the median marker expression per sample
expr_median_sample_tbl <- data.frame(sample_id = metadata_sc$sample_id, expr) %>%
  group_by(sample_id) %>% summarize_all(list(median)) %>% as.data.frame()
rownames(expr_median_sample_tbl) <- expr_median_sample_tbl$sample_id
expr_median_sample_tbl <- expr_median_sample_tbl[,-1]

#compute the PCA
drPCA <- dim.reduction(expr_data = expr_median_sample_tbl, metadata = md, reduction_method = "PCA")
dr.plotting(drPCA$dr_melted, dr_calculated = "PCA", color_by = "condition", size = 3, output_type = "png", labels = "patient_id")

#prepare data and annotations for heatmap
annotation_col <- data.frame(row.names = rownames(expr_median_sample_tbl), condition = md$condition)
color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100) #colors for the heatmap (expression values)
x <- c()
set.seed(3); for(i in unique(md$condition)) x[i] <- sample(colors_palette, 1)
annotation_colors <- list(condition = x) #change the name (ie "condition") for according to annotation in heatmap

pheatmap(t(expr_median_sample_tbl), color = color, display_numbers = FALSE,
         number_color = "black", fontsize_number = 5, annotation_col = annotation_col,
         annotation_colors = annotation_colors, clustering_method = "average")

## PCA and heatmap on single cell expression
sub_idx_heat <- sub.samples.idx(data = mdsc, colname_samples = "sample_id", samples_names = md$sample_id, 
                         subsampling = 1000, set.seed = 1234) #select how many cells to downsample per-sample
heat_expr<- mdsc[sub_idx_heat,]

#compute PCA (single cell)
ggdfPCA_1000 <- dim.reduction(expr_data = heat_expr[,4:ncol(heat_expr)], metadata = heat_expr[,1:3], 
                              reduction_method = "PCA")

#PCA colored by condition or patient -> batch effect?
dr.plotting(ggdfPCA_1000$dr_melted, dr_calculated = "PCA", color_by = "condition", output_type = NULL)
dr.plotting(ggdfPCA_1000$dr_melted, dr_calculated = "PCA", color_by = "patient_id", output_type = NULL)

#heatmap (single cell)
pheatmap(t(heat_expr[,4:ncol(heat_expr)]), color = color, display_numbers = FALSE,
         number_color = "black", fontsize_number = 5, annotation_col = heat_expr[,2:3], show_colnames = F, 
         annotation_colors = annotation_colors, clustering_method = "average", treeheight_row = 0, treeheight_col = 0)


#################################################
### Cell clustering #############################
#################################################
#### FlowSOM: SOM and MST ####
fsom <- fsom.clustering(fcs, markers_to_use = "surface_markers", markers_to_plot = "tree", set.seed = 1234)

#### FlowSOM: metaclustering ####
metaclusters <- fsom.metaclustering(fsom = fsom, num_clusters_metaclustering = 40, plotting = T, set.seed = 1234)

#heatmap
cluster_heatmap(expr = expr[, surface_markers], expr_saturated = expr01[, surface_markers],
                                cell_clusters = metaclusters$metaclusters)

#### Dimensional reduction: PCA, tSNE and UMAP ####
#combine single-cell data with SOM clusters
mdsc_som <- data.frame(mdsc[,1:3], SOM = as.factor(metaclusters$metaclusters), mdsc[,4:ncol(mdsc)])

#subsampling
sub_idx_som <- sub.samples.idx(data = mdsc_som, colname_samples = "sample_id", samples_names = md$sample_id, subsampling = 1000, set.seed = 1234)
sel_expr <- mdsc_som[sub_idx_som,]

## Calculate dimensional reductions (simultaneously)
dr <- dim.reduction(expr_data = sel_expr[,5:ncol(sel_expr)], metadata = sel_expr[,1:4], reduction_method = "all", set.seed = 1234)

#plot PCA, tSNE or UMAP -> change "dr_calculated"
dr.plotting(dr$dr_melted, dr_calculated = "tSNE", output_type = "png", color_by = "SOM")
dr.plotting(dr$dr, dr_calculated = "UMAP", output_type = "png", color_by = "CD62L")

for(facet in c("patient_id", "condition", "sample_id")){
  dr.plotting(dr$dr_melted, dr_calculated = "PCA", output_type = "png", color_by = "SOM", facet_by = facet, 
              output_name = paste0(facet, "_"))
}

#################################################
### FCS exporting and cluster analysis ##########
#################################################
#### Generate FCS files with dimensional reduction data ####
to_export <- data.frame(sample_id = as.numeric(sel_expr$sample_id), patient_id = as.numeric(sel_expr$patient_id), 
                        condition = as.numeric(sel_expr$condition), cluster_som = as.numeric(sel_expr$SOM),
                        expr_no_transf[sub_idx_som,], dr$dr[,grepl("PCA|tSNE|UMAP", colnames(dr$dr))])

outFile <- file.path(getwd(), "alltubeTh.fcs")
write.FCS(as_flowFrame(as.matrix(to_export), source.frame = NULL), outFile)

# ++++++++++++++++++++++++++++
# external step >>> manual analysis in a flow cytometry software to identify clusters >>> Excel file
# ++++++++++++++++++++++++++++

#### Merge clusters according external analysis ####
# cluster_merging1_filename <- "/mnt/beegfs/cdalisvam/paper/Th_new_clustering.xlsx"
cluster_merging1 <- read_excel("../../Th_new_clustering.xlsx", col_types = "text")
head(cluster_merging1)

#add new cluster names to the original metadata matrix 
cell_clustering1m <- FindReplace(data.frame(cell_clustering1m = as.factor(metaclusters$metaclusters)), Var = "cell_clustering1m", 
                                 replaceData = cluster_merging1, 
                                 from = "original_cluster", to = "new_cluster", vector = T) %>% as.factor()

cell_clustering1m_plotStars <- FindReplace(data.frame(cell_clustering1m = as.factor(metaclusters$plotStars_value)), Var = "cell_clustering1m", 
                                 replaceData = cluster_merging1, 
                                 from = "original_cluster", to = "new_cluster", vector = T) %>% as.factor()

dr$dr_melted$cell_clustering1m <- FindReplace(data.frame(cell_clustering1m = as.factor(dr$dr$SOM)), 
                                    Var = "cell_clustering1m", replaceData = cluster_merging1, 
                                    from = "original_cluster", to = "new_cluster", vector = T) %>% as.factor()

#draw PCA, tSNE, MST and heatmap colored by new merged clusters
dr.plotting(dr$dr_melted, dr_calculated = "UMAP", color_by = "cell_clustering1m", output_type = NULL)

PlotStars(fsom, backgroundValues = cell_clustering1m_plotStars,  
  backgroundColor = alpha(colors_palette, alpha = 0.4))

cluster_heatmap(expr = expr[, surface_markers],
                                expr_saturated = expr01[, surface_markers], cell_clusters = cell_clustering1m)


#################################################
### Subclustering ###############################
#################################################
#### Generation of a new matrix for subclustering ####
metadata_sc <- cbind(mdsc, FlowSOM = cell_clustering1m)
head(metadata_sc)

#remove FSC_H and SSC_H that will not be used for further analysis
matrix_clusters <- as.matrix(metadata_sc[metadata_sc$FlowSOM == "lymphocytes", surface_markers])
head(matrix_clusters)

rngL <- colQuantiles(matrix_clusters, probs = c(0.01, 0.99))
matrix_clusters01 <- t((t(matrix_clusters) - rngL[, 1]) / (rngL[, 2] - rngL[, 1]))
matrix_clusters01[matrix_clusters01 < 0] <- 0
matrix_clusters01[matrix_clusters01 > 1] <- 1

#### Cell clustering and dimensional reduction ####
#select population for subclustering
metadata_scL <- metadata_sc[metadata_sc$FlowSOM == "lymphocytes",]
rownames(metadata_scL) <- NULL #reset rownames after cell selection (for later subsampling)

expr_no_transfL <- expr_no_transf[metadata_sc$FlowSOM == "lymphocytes",]

#lymphocytes subclustering
cell.count.bx(metadata_scL, counts_by = "sample_id", metadata = md)

## PCA reduction (median values)
#get the median marker expression per sample
expr_median_sampleL <- data.frame(sample_id = metadata_scL$sample_id, metadata_scL[,surface_markers]) %>%
  group_by(sample_id) %>% summarize_all(funs(median)) %>% as.data.frame()
rownames(expr_median_sampleL) <- expr_median_sampleL$sample_id
expr_median_sampleL <- expr_median_sampleL[,-1]

#compute PCA and visualization
PCA_outL <- dim.reduction(expr_data = expr_median_sampleL, metadata = md, reduction_method = "PCA")
dr.plotting(PCA_outL$dr_melted, dr_calculated = "PCA", color_by = "condition", output_type = "png", size = 6)

## FlowSOM
fsomL <- fsom.clustering(data = metadata_scL[,surface_markers], markers_to_use = "surface_markers", set.seed = 1234)
metaclustersL <- fsom.metaclustering(fsom = fsomL, num_clusters_metaclustering = 40, plotting = T, set.seed = 1234)

cluster_heatmap(expr = metadata_scL[,surface_markers], expr_saturated = matrix_clusters01[,surface_markers], 
                                cell_clusters = metaclustersL$metaclusters)

metadata_scL$FlowSOM_L <- metaclustersL$metaclusters

## PCA and tSNE
#subsampling
sub_idxL <- sub.samples.idx(metadata_scL, colname_samples = "sample_id", samples_names = md$sample_id, 
                            subsampling = 1000, set.seed = 1234)
sel_exprL <- metadata_scL[sub_idxL,]

#dimensional reduction
drL <- dim.reduction(sel_exprL[,surface_markers], metadata = sel_exprL[,!(colnames(sel_exprL) %in% surface_markers)], 
                     reduction_method = "all", set.seed = 1234)

#change "dr_calculated" to plot each dimensional reduction
gg_list <- list() #list to store plots
for(marker in surface_markers) gg_list[[marker]] <- dr.plotting(drL$dr, dr_calculated = "PCA", color_by = marker, output_type = NULL)
do.call("grid.arrange", c(gg_list, nrow = 2)) #show all plots

dr.plotting(drL$dr_melted, dr_calculated = "tSNE", color_by = "FlowSOM_L", output_type = NULL, facet_by = "patient_id")

## Generate new FCS with dimensional reduction data for lymphocytes
to_export <- data.frame(sample_id = as.numeric(sel_exprL$sample_id), patient_id = as.numeric(sel_exprL$patient_id), 
                        condition = as.numeric(sel_exprL$condition), cluster_somL = as.numeric(sel_exprL$FlowSOM_L),
                        expr_no_transfL[sub_idxL,], drL$dr[,grepl("PCA|tSNE|UMAP", colnames(drL$dr))])

outFileL <- file.path(getwd(), "alltubeThclustL.fcs")
write.FCS(as_flowFrame(as.matrix(to_export), source.frame = NULL), outFileL)

# ++++++++++++++++++++++++++++
# external step >>> manual analysis in a flow cytometry software to identify clusters >>> Excel file
# ++++++++++++++++++++++++++++

#### Identification and representation of multiple lymphocyte subsets ####
## Merge clusters according external analysis
cluster_merging1L <- read_excel("../../../Th_new_clusteringL.xlsx", col_types = "text")
head(cluster_merging1L)

cell_clustering1mL <- FindReplace(data.frame(cell_clustering1m = as.factor(metaclustersL$metaclusters)), Var = "cell_clustering1m", 
                                 replaceData = cluster_merging1L, 
                                 from = "original_cluster", to = "new_cluster", vector = T) %>% as.factor()

cell_clustering1mL_plotStars <- FindReplace(data.frame(cell_clustering1m = as.factor(metaclustersL$plotStars_value)), Var = "cell_clustering1m", 
                                  replaceData = cluster_merging1L, 
                                  from = "original_cluster", to = "new_cluster", vector = T) %>% as.factor()

drL$dr_melted$cell_clustering1mL <- FindReplace(data.frame(cell_clustering1m = as.factor(drL$dr$FlowSOM_L)), 
                                              Var = "cell_clustering1m", replaceData = cluster_merging1L, 
                                              from = "original_cluster", to = "new_cluster", vector = T) %>% as.factor()


#FlowSOM
PlotStars(fsomL, backgroundValues = cell_clustering1mL_plotStars,  backgroundColor = alpha(colors_palette, alpha = 0.5))

#Dimensional reduction plotting -> change to PCA, tSNE or UMAP
dr.plotting(drL$dr_melted, dr_calculated = "tSNE", color_by = "cell_clustering1mL", output_type = NULL, facet_by = "patient_id")

cluster_heatmap(expr = matrix_clusters[,surface_markers], 
                                expr_saturated = matrix_clusters01[,surface_markers],
                                cell_clusters = cell_clustering1mL)

## Circular hyerarchical clustering tree
circ.tree.selectNodes(exprs = matrix_clusters, exprs_saturated = matrix_clusters01, cell_clusters = cell_clustering1mL)
circ.tree(exprs = matrix_clusters, exprs_saturated = matrix_clusters01, cell_clusters = cell_clustering1mL, 
          dendro_labels = F, nodes = c(29,33,35))


## Heatmap for single cells
#subsample samples
sub_idxL_heat <- sub.samples.idx(metadata_scL, colname_samples = "sample_id", samples_names = md$sample_id, 
                            subsampling = 1000, set.seed = 1234)
sel_exprL_heat <- metadata_scL[sub_idxL_heat,]

#prepare heatmap
annotation_col <- data.frame(row.names = rownames(sel_exprL_heat), sel_exprL_heat[,!(colnames(sel_exprL_heat) %in% surface_markers)])
annotation_colors <- list(condition = c(BM = colors_palette[1], PB = colors_palette[2]))

pheatmap(t(sel_exprL_heat[,surface_markers]), annotation_col = annotation_col, clustering_method = "average",
         color = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100), display_numbers = FALSE,
         number_color = "black", fontsize_number = 5, annotation_colors = annotation_colors, 
         show_colnames = F, treeheight_col = 0, treeheight_row = 0)


#################################################
### Data export and statistical analysis ########
#################################################
#### Visualization of cellular abundance and exporting results in a tabular format ####
## Visualization
prop_tableL <- barplot.cell.pops(cell_clustering1mL, metadata = metadata_scL, colname_sampleID = "sample_id")

## Data exporting
dataset <- data.frame(md, as.data.frame.matrix(t(prop_tableL)))
write.table(dataset, file = "resultsL.txt", sep = "\t")

# library(haven) #format for SPSS
# write_sav(dataset, "resultsL.sav")

#### Evaluation of normal distribution and statistical analysis ####
## Normality tests
#Shapiro test (p<0.05, distribution not normal)
lshap <- lapply(dataset[,5:ncol(dataset)], shapiro.test) #globally
lshap[[1]]

lres <- sapply(lshap, `[`, c("statistic","p.value")) #by cell type
t(lres)

#QQ plots
qqnorm(dataset$CD8_naive_CCR4p) #not normal distribution
qqnorm(dataset$CD4_tm_DP) #normal distribution

#densityplot
ggdensity(dataset, x = "CD8_EM_CD27p",
          add = "mean", rug = TRUE,
          color = "condition", fill = "condition",
          palette = colors_palette) + xlim(-10,70)

## Correlation between condittions
cor.plot.conditions(dataset, colname_condition = "condition", colname_patientID = "patient_id", 
                    conditions = c("SP", "MO"), metadata = md)


## Boxplot visualization
ggdf_L_calcula <- melt(dataset, id.vars = c("file_name", "sample_id", "condition", "patient_id"), 
                       variable.name = "cluster", value.name = "proportion")

resultskwSP_BM <- col_kruskalwallis(dataset[,5:ncol(dataset)], dataset$condition)
KWsig <- rownames(resultskwSP_BM[resultskwSP_BM$pvalue<0.1,])

#those significant populations
ggboxplot(ggdf_L_calcula[ggdf_L_calcula$cluster %in% KWsig,], x = "condition", y = "proportion",
          color = "condition", fill = "condition", alpha = 0.5,
          add = "jitter", facet.by = "cluster", short.panel.labs = FALSE, nrow = 2)+ 
  geom_point(aes(x = condition, y = proportion, color = condition,
                 shape = patient_id), alpha = 0.8, position = position_jitterdodge()) +
  scale_y_continuous(trans = "log10") + #log transformation, evaluate differences visually 
  # theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), strip.text = element_text(size = 10))+
  scale_color_manual(values = colors_palette) + scale_fill_manual(values = colors_palette) +
  scale_shape_manual(values = 1:length(unique(ggdf_L_calcula$patient_id)))+
  stat_compare_means(aes(group = condition), method = "kruskal.test", label = "p.format", label.y = 2)

#all populations
ggboxplot(ggdf_L_calcula, x = "condition", y = "proportion",
             color = "condition",fill = "condition",alpha = 0.5, 
             add = "jitter", facet.by = "cluster", short.panel.labs = FALSE, nrow = 3)+ 
  geom_point(aes(x = condition, y = proportion, color = condition,
                 shape = patient_id), alpha = 0.8, position = position_jitterdodge()) +
  scale_y_continuous(trans = "log10") +
  # theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), strip.text = element_text(size = 10))+
  scale_color_manual(values = colors_palette) + scale_fill_manual(values = colors_palette) +
  scale_shape_manual(values = 1:length(unique(ggdf_L_calcula$patient_id)))+
  stat_compare_means(aes(group = condition), method = "kruskal.test", label = "p.format", label.y = 2 )

#prepare final table
final_metadata <- data.frame()
subclust_samples <- list(lymphocytes = cell_clustering1mL) #those cellular populations used for previous subclustering 

for(population in levels(metadata_sc$FlowSOM)){
  aux <- metadata_sc[metadata_sc$FlowSOM == population,]
  
  if(population %in% names(subclust_samples)){
    aux$final_cluster <- as.character(subclust_samples[[population]])
  }else{
    aux$final_cluster <- as.character(aux$FlowSOM)
  }
  
  final_metadata <- rbind(final_metadata, aux)
}
final_metadata$final_cluster <- as.factor(final_metadata$final_cluster)

#apply on dimensional reduction
mm <- match(rownames(dr$dr), rownames(final_metadata))
dr$dr$cell_clustering_final <- final_metadata[mm, "final_cluster"]

dr.plotting(dr$dr, dr_calculated = "PCA", color_by = "cell_clustering_final", facet_by = "patient_id", output_type = NULL)

#export final FCS file with definitive clusters
to_export <- data.frame(sample_id = as.numeric(final_metadata$sample_id)[sub_idx_som], 
                        patient_id = as.numeric(final_metadata$patient_id)[sub_idx_som], 
                        condition = as.numeric(final_metadata$condition)[sub_idx_som], 
                        final_cluster = as.numeric(final_metadata$final_cluster)[sub_idx_som], 
                        expr_no_transf[sub_idx_som,], dr$dr[,grepl("PCA|tSNE|UMAP", colnames(dr$dr))])

outFilefinal<- file.path(getwd(), "Thclustfinal.fcs")
write.FCS(as_flowFrame(as.matrix(to_export), source.frame = NULL), outFilefinal)

#differential cell population abundance
props_final <- barplot.cell.pops(cell_clusters = final_metadata$final_cluster, metadata = final_metadata[,1:3], colname_sampleID = "sample_id")

#export final results for external analysis
to_export2 <- cbind(md, as.data.frame.matrix(t(props_final)))
write.table(to_export2, file="results_final.txt",sep="\t")

#test again for normal distribution
dataset_final <- to_export2

lshap_final <- lapply(dataset[-(1:4)], shapiro.test) #Shapiro test
lshap_final[[1]]

lres_final <- sapply(lshap_final, `[`, c("statistic","p.value"))
t(lres_final)

qqnorm(dataset_final$granulocytes) #QQplot

ggdensity(dataset_final, x = "erythroblasts",#densityplot
          add = "mean", rug = TRUE,
          color = "condition", fill = "condition",
          palette = colors_palette) + xlim(-10,70)

#correlation between cell distributions BM/PB
cor.plot.conditions(data = dataset_final, colname_condition = "condition", colname_patientID = "patient_id", 
                    conditions = c("PB", "BM"), metadata = md)

#comparisons between conditions
props_final_calcula <- melt(dataset_final, id.vars = c("file_name", "sample_id", "condition", "patient_id"), 
                       variable.name = "cluster", value.name = "proportion")

resultskwSP_BM_final <- col_kruskalwallis(dataset_final[,5:ncol(dataset_final)], dataset_final$condition)
KWsig_final<-rownames(resultskwSP_BM_final[resultskwSP_BM_final$pvalue<0.05,])

ggboxplot(props_final_calcula[props_final_calcula$cluster %in% KWsig_final,], x = "condition", y = "proportion",
          color = "condition", 
          fill = "condition", alpha = 0.5, 
          add = "jitter", facet.by = "cluster", short.panel.labs = FALSE, nrow = 2)+ 
  geom_point(aes(x = condition, y = proportion, color = condition,
                 shape = patient_id), alpha = 0.8, position = position_jitterdodge()) +
  scale_y_continuous(trans = "log10") + #log transformation, evaluate differences visually 
  # theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), strip.text = element_text(size = 10))+
  scale_color_manual(values = colors_palette) + scale_fill_manual(values = colors_palette) +
  scale_shape_manual(values = c(1:length(props_final_calcula$patient_id)))+
  stat_compare_means(aes(group = condition), method = "kruskal.test", label = "p.format", label.y = 2)
```
