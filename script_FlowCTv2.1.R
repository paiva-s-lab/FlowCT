## INFORMATION ##########################################################################
# FlowCT version ---> 2.2
# backbone script version ---> 1.1
##########################################################################################

### Environment setting ##################################################################
## Load packages and functions
# Sys.setenv(http_proxy  = "http://proxy.unav.es:8080")
# Sys.setenv(https_proxy = "http://proxy.unav.es:8080")
# devtools::install_github("jgarces02/FlowCT@devel", auth_token = "21ea9880f944d42755479e54a5b19ddd00fe17f6")
library(FlowCT.v2)

## Working directory
setwd("G:/Mi unidad/Proyectos/FlowCT_Ciro/FlowCT.v2/results/")
# sapply(grep(list.files(path = "../src/", full.names = T), pattern = 'script|README|loading|old', inv = T, value = T), source)

### unify all FCS nomenclatures ##########################################################
# unify.FCSheaders(directory = "../data/", pattern = "fcs", fix = F) #bug: update data.table::melt ---> done
# unify.FCSheaders(directory = "../data/", pattern = "fcs", fix = T, select.freq = 2)

### Prepare metadata and fcs.SCE object ##################################################
(filenames <- list.files(pattern = "fcs", path = "../data/"))

md <- data.frame(filename = filenames, #mandatory
                 sample_id = 1:length(filenames),
                 condition = sapply(filenames,function(x) strsplit(x, split = "_|\\.")[[1]][1]),
                 patient_id = sapply(filenames,function(x) strsplit(x, split = "_|\\.")[[1]][2]))
(md)

## FCS reading and transforming
system.time(fcs <- fcs.SCE(directory = "../data/", pattern = "fcs", events = 1000, metadata = md,
                           transf.cofactor = 500, project.name = "paired"))

## adjust maker names
marker.names(fcs)
new_names <- c("FSC_A", "FSC_H", "SSC_A", "SSC_H", "CD62L", "CXCR3", "CD8", "CCR4", "CCR6", "CD4", "CD45", "CD27")
fcs <- marker.names(fcs, new.names = new_names)

## select markers for downstrean analysis
surface_markers <- c("CD62L", "CXCR3", "CD8", "CCR4", "CCR6", "CD4", "CD45", "CD27")
physical_markers <- c("FSC_A", "SSC_A")


### QC and doublets removal ##############################################################
fcs <- qc.and.removeDoublets(fcs, physical.markers = c("FSC_A", "FSC_H", "SSC_A", "SSC_H"),
                             return.fcs = F)


### Normalization and data alignment #####################################################
multidensity(fcs.SCE = fcs, assay.i = "transformed", subsampling = 100)

fcs <- normalization.flw(fcs.SCE = fcs, marker.to.norm = c("CCR6", "CCR4"),
                         norm.method = "warp")

multidensity(fcs.SCE = fcs, assay.i = "normalized", subsampling = 1000)
multidensity(fcs, assay.i = 2, ridgeline.lim = 10, color.by = "filename",
             show.markers = c("CD62L", "CD4"), interactive = T)


### Descriptive and exploratory analysis #################################################
## boxplot with cell numbers by condition
cell.count.bx(fcs.SCE = fcs, x.axis = "condition") #bug: add option to chose assay ---> done

## dimensinal reduction and heatmap with median values
median.dr(fcs, color.by = "filename")
median.heatmap(fcs, not.metadata = c("sample_id", "filename"))

## PCA and heatmap on single cell expression
fcs100 <- sub.samples(fcs.SCE = fcs, subsampling = 100)

fcs100 <- dim.reduction(fcs100, dr.method = c("tSNE", "pca"))
dr.plotting(fcs100, plot.dr = "tsne", color.by = "condition") #bug: install scattermore ---> done

sc.heatmap(fcs100, subsampling = 100, markers.to.use = surface_markers) #bug: Floating point exception (???) ---> done!


### Clustering ###########################################################################
## FlowSOM
fsom <- fsom.clustering(fcs.SCE = fcs, markers.to.use = surface_markers, k.metaclustering = 40,
                        markers.to.plot = "tree_metaclustering") #bug: it's calculating clustering with a fcs.SCE incorrectly generated (w/o "normalized" assay) > why?? It must bring error!
fcs$SOM <- fsom$metaclusters

median.heatmap(fcs.SCE = fcs, assay.i = "normalized", cell.clusters = fcs$SOM)

## dimensional reduction
fcs1000 <- sub.samples(fcs.SCE = fcs, subsampling = 1000)

fcs1000 <- dim.reduction(fcs1000, dr.method = c("PCA", "UMAP", "tsNe"))

dr.plotting(fcs1000, plot.dr = "tSNE", color.by = "condition")
dr.plotting(fcs1000, plot.dr = "UMAP", color.by = "patient_id")
dr.plotting(fcs1000, plot.dr = "PCA", color.by = "SOM", facet.by = "condition")

## FCS exporting and cluster analysis
## Generate FCS files with dr data
export.metaFCS(fcs.SCE = fcs1000, output.name = "hola")
export.metaFCS(fcs.SCE = fcs1000, separate.fcs = T, output.suffix = "s")

# ++++++++++++++++++++++++++++
# external step >>> manual analysis in a flow cytometry software to identify clusters >>> Excel file
# ++++++++++++++++++++++++++++

## rename and merge clusters according external analysis
replacedata <- readxl::read_excel("../data/Th_new_clustering.xlsx", col_types = "text")

fcs$SOM_named <- clusters.rename(fcs$SOM, cluster = replacedata$original_cluster, name = replacedata$new_cluster)
fcs1000$SOM_named <- clusters.rename(fcs1000$SOM, cluster = replacedata$original_cluster, name = replacedata$new_cluster)
fsom$plotStars_value_named <- clusters.rename(fsom$plotStars_value, cluster = replacedata$original_cluster, name = replacedata$new_cluster)

## draw PCA, tSNE, MST and heatmap colored by new merged clusters
dr.plotting(fcs1000, plot.dr = "UMAP", color.by = "SOM_named")

PlotStars(fsom$fsom, backgroundValues = fsom$plotStars_value_named,
          backgroundColor = alpha(div.colors(40), alpha = 0.4), starBg = NULL)

median.heatmap(fcs.SCE = fcs1000, assay.i = "normalized", cell.clusters = fcs1000$SOM_named)


### Subclustering ########################################################################
## prepare fcs.SCE object
fcsL <- fcs[,fcs$SOM_named == "lymphocytes"]

## exploratory analysis
cell.count.bx(fcsL, assay.i = "normalized", x.axis = "condition")

median.heatmap(fcsL, not.metadata = c("sample_id", "filename"))

## FlowSOM
fsomL <- fsom.clustering(fcs.SCE = fcsL, markers.to.use = surface_markers, k.metaclustering = 40)
fcsL$SOM_L <- fsomL$metaclusters

median.heatmap(fcs.SCE = fcsL, assay.i = "normalized", cell.clusters = fcsL$SOM_L)

## DR
fcsL <- dim.reduction(fcsL, dr.method = c("tsne", "pca", "umap"), markers.to.use = surface_markers)

dr.plotting(fcsL, plot.dr = "Umap")
dr.plotting(fcsL, plot.dr = "tSNE", color.by = "SOM_L")

## Generate FCS files with subclustering data
export.metaFCS(fcs.SCE = fcsL, output.name = "subclust1.vSCE.fcs")

# ++++++++++++++++++++++++++++
# external step >>> manual analysis in a flow cytometry software to identify clusters >>> Excel file
# ++++++++++++++++++++++++++++

## rename and merge clusters according external analysis
replacedataL <- readxl::read_excel("../data/Th_new_clusteringL.xlsx", col_types = "text")

fcsL$SOM_L_named <- clusters.rename(fcsL$SOM_L, cluster = replacedataL$original_cluster, name = replacedataL$new_cluster)
fsomL$plotStars_value_named <- clusters.rename(fsomL$plotStars_value, cluster = replacedataL$original_cluster, name = replacedataL$new_cluster)

## draw PCA, tSNE, MST and heatmap colored by new merged clusters
dr.plotting(fcsL, plot.dr = "tSNE", color.by = "SOM_L_named")

PlotStars(fsomL$fsom, backgroundValues = fsomL$plotStars_value_named,
          backgroundColor = alpha(div.colors(40), alpha = 0.4), starBg = NULL)

median.heatmap(fcs.SCE = fcsL, assay.i = "normalized", cell.clusters = fcsL$SOM_L_named)

## Circular hyerarchical clustering tree
circ.tree(fcs.SCE = fcsL, cell.clusters = fcsL$SOM_L_named, nodes = "display")
circ.tree(fcs.SCE = fcsL, cell.clusters = fcsL$SOM_L_named, nodes = c(29, 37))

## heatmap for single cells
sc.heatmap(fcsL, assay.i = "normalized", not.metadata = c("cell_ID", "filename", "SOM", "SOM_named", "SOM_L"))


### cluster identif (misc.) ##############################################################
## differential dotplots
dotplot.DE(fcs.SCE = fcsL, markers.to.use = surface_markers, clusters.named = "SOM_L_named")


### data exporting #######################################################################
## delete non-interesting populations
delete_pops <- "unclassified"
fcs1000_rm <- remove.pop(fcs1000, clusters.named = "SOM_named", population = delete_pops)
fcsL_rm <- remove.pop(fcsL, clusters.named = "SOM_L_named", population = delete_pops)

## combine clustering and subclustering
fcs_final <- combine.subclusterings(initial.fcs.SCE = fcs1000, clusters.named = "SOM_named",
                                    subclustering.fcs.SCE = list(fcsL_rm))

## final FCS files with named clusters and DR info
export.metaFCS(fcs.SCE = fcs_final, output.name = "prueba_final4C.fcs")

## cellular proportions (or raw counts)
prop_tableL <- barplot.cell.pops(fcs.SCE = fcs_final, cell.clusters = fcs_final$SOM_named_final,
                                 count.by = "sample_id", facet.by = "condition",
                                 return.mode = "percentage")

dataset <- data.frame(md, as.data.frame.matrix(t(prop_tableL)))
write.table(dataset, file = "results.txt", sep = "\t")

## median values (from non-transformed data)
med <- median.values(fcs.SCE = fcs_final)
# med <- median.values(fcs.SCE = fcs.SCE_final, var = "SOM_named_final")
dataset2 <- data.frame(md, med)
write.table(dataset2, file = "results_median.txt", sep = "\t")

### statistics ###########################################################################
## Correlation between (two) conditions
corplot.conditions(fcs.SCE = fcsL, cell.clusters = fcsL$SOM_L_named,
                   condition.column = "condition")

corplot.conditions(fcs.SCE = fcs_final, cell.clusters = fcs_final$SOM_named_final,
                   condition.column = "condition")

## boxplot
boxplot.cell.clustering(fcs.SCE = fcsL, cell.clusters = fcsL$SOM_L_named, facet = T,
                        return.stats = F, plot.only.sig = c(T, 0.1), facet.free.scale = "free")
sig <- boxplot.cell.clustering(fcs.SCE = fcs_final, cell.clusters = fcs_final$SOM_named_final,
                               return.stats = T)

## Dumbbell plot
dumbPlot.cell.clustering(fcs.SCE = fcs_final, cell.clusters = fcs_final$SOM_named_final,
                         return.stats = F, condition.column = "condition")

## diffdots plot
diffdots.cell.clustering(fcs.SCE = fcs_final, cell.clusters = fcs_final$SOM_named_final,
                         return.stats = F, condition.column = "condition")
