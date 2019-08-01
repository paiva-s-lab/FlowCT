#' plot_clustering_heatmap_wrapper
#'
#' This function draws a heatmap collapsing identified clusters and showing cell percentaje for each one.
#' @param expr A data.frame with marker expression values for each cell
#' @param expr_saturated A data.frame with marker expression values for each cell where those expression values greater or smaller that 1 and 0, respectively, must be coherced to 1 and 0.
#' @param cell_clusters A vector with clusters identified in a previous clustering method.
#' @keywords heatmap
#' @keywords clusters
#' @keywords cell percentajes
#' @export
#' @import dplyr
#' @import stats
#' @importFrom ggplot2 ggsave
#' @importFrom grDevices colorRampPalette
#' @importFrom pheatmap pheatmap
#' @examples
#' \dontrun{
#' plot_clustering_heatmap_wrapper(expr = expr[,surface_markers],
#'     expr_saturated = expr01[,surface_markers], cell_clusters = df$FlowSOM)}

cluster.heatmap <- function(expr, expr_saturated, cell_clusters, cluster_merging = NULL){
  #calculate the median expression
  expr_median <- data.frame(expr, cell_clusters = cell_clusters) %>%
    group_by(cell_clusters) %>% summarize_all(list(median))
  
  expr_saturated_median <- data.frame(expr_saturated, cell_clusters = cell_clusters) %>%
    group_by(cell_clusters) %>% summarize_all(list(median))
  
  #calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clusters))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  
  #sort cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  expr_heat <- as.matrix(expr_saturated_median[, colnames(expr_saturated)])
  rownames(expr_heat) <- expr_saturated_median$cell_clusters
  
  #colors for the heatmap
  color_heat <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
  legend_breaks <- seq(from = 0, to = 1, by = 0.1)
  labels_row <- paste0(expr_saturated_median$cell_clusters, " (", clustering_prop ,"%)")
  
  #annotation for original clusters
  annotation_row <- data.frame(Cluster = factor(expr_saturated_median$cell_clusters))
  rownames(annotation_row) <- rownames(expr_heat)
  colors_palette1 <- colors_palette[1:nlevels(annotation_row$Cluster)]
  names(colors_palette1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = colors_palette1)
  
  #annotation for merged clusters
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Cluster_merging <- cluster_merging$new_cluster
    colors_palette2 <- colors_palette[1:nlevels(cluster_merging$new_cluster)]
    names(colors_palette2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Cluster_merging <- colors_palette2
  }
  
  p <- pheatmap(expr_heat, color = color_heat, cluster_cols = FALSE,
           cluster_rows = cluster_rows, labels_row = labels_row,
           display_numbers = FALSE, number_color = "black",
           fontsize = 8, fontsize_number = 6, legend_breaks = legend_breaks,
           annotation_row = annotation_row, annotation_colors = annotation_colors)

  ## check X11 active to redirect output
  if("try-error" %in% class(suppressWarnings(try(x11(), silent = T)))){
    cat("X11 is not active, percentaje heatmap is saved as -> pctHeatmap.", deparse(substitute(expr)), ".jpg\n", sep = "")
    suppressMessages(ggsave(paste0("pctHeatmap.", gsub("\\[.*.", "", deparse(substitute(expr))), ".jpg"), device = "jpeg", plot = p))
  }else{
    p
  }
}


#' dr.plotting
#'
#' This function plots a selected dimensional reduction from those previously calculated with \code{\link[FlowCT:dim.reduction]{FlowCT::dim.reduction()}}.
#' @param data A data.frame generated with \code{\link[FlowCT:dim.reduction]{FlowCT::dim.reduction()}}.
#' @param dr_calculated Type of dimensional reduction for plotting. Possible values are "PCA", "tSNE" or "UMAP".
#' @param color_by Variable for coloring the plot. Default = "expression"
#' @param facet_by Variable for faceting the plot. Default = \code{NULL}
#' @param labels Variable for labeling points. Default = \code{NULL}
#' @param size Point size. Default = 0.5
#' @param output_name Preffix to add to the final saved file. Default = ""
#' @param output_type Format for saving plot. Possible values are \code{NULL} (Default, i.e. the plot is shown in the \emph{Plots panel}), "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf".
#' @keywords dimensional reduction
#' @keywords tSNE
#' @keywords PCA
#' @keywords UMAP
#' @export
#' @import ggplot2
#' @examples
#' \dontrun{
#' dr.plotting(dr$dr_melted, dr_calculated = "PCA", color_by = "patient",
#'     output_type = NULL)
#' dr.plotting(dr$dr_melted, dr_calculated = "tSNE", color_by = "patient",
#'     facet_by = "condition", output_type = "png")}

dr.plotting <- function(data, dr_calculated, color_by = "expression", facet_by = NULL, 
                        labels = NULL, size = 0.5, output_name = "", output_type = NULL){
  
  g <- ggplot(data, aes(x = eval(parse(text = paste0(dr_calculated, "1"))), y = eval(parse(text = paste0(dr_calculated, "2"))))) +
    geom_point(aes(color = eval(parse(text = color_by))), size = size) + 
    xlab(paste0(dr_calculated, 1)) + ylab(paste0(dr_calculated, 2)) +
    ggtitle(paste0(dr_calculated, " reduction, colored by ", color_by)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA))
  
  if(is.factor(data[,color_by])){
    g <- g + scale_color_manual(values = colors_palette, name = color_by) +
      guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
  }else{
    g <- g + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50), name = color_by)
  }
  
  if(!is.null(facet_by)){
    g <- g + facet_wrap(~ eval(parse(text = facet_by)))
  }
  
  if(!is.null(labels)){
    g <- g + geom_text(aes(label = eval(parse(text = labels))), nudge_y = 0.05)
  }
  
  if(is.null(output_type)){
    print(g)    
  }else{
    suppressMessages(ggsave(plot = g, filename = paste0(output_name, "dr.", dr_calculated, "_col.", color_by, ".", output_type), device = output_type))
  }
  
  cat(paste0(dr_calculated, " reduction colored by ", color_by, 
             ifelse(is.null(facet_by), "", paste0(" and faceted by ", facet_by)), 
             ifelse(is.null(output_type), "", paste0(". Saved as -> ", paste0(output_name, "dr.", dr_calculated, "_col.", color_by, ".", output_type)))), "\n")
  return(g)
}


#' circ.tree.selectNodes
#'
#' This function plots a circular dendrogram for clusters previously calculated in order to select those nodes for coloring in later \code{\link[FlowCT:circ.tree]{FlowCT::circ.tree()}}.
#' For additional information go to \href{https://guangchuangyu.github.io/software/ggtree/documentation/}{\code{\link{ggtree}}} package.
#' @param exprs A data.frame with marker expression values for each cell
#' @param exprs_saturated A data.frame with marker expression values for each cell where those expression values greater or smaller that 1 and 0, respectively, must be coherced to 1 and 0.
#' @param cell_clusters A vector with clusters identified in a previous clustering method.
#' @param dist_method The distance measure to be used. Possible values are "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param hclust_method The agglomeration method to be used. Possible values are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".
#' @keywords circular tree
#' @keywords dendrogram
#' @keywords node
#' @export
#' @import dplyr
#' @import stats
#' @import ggtree
#' @importFrom ape as.phylo
#' @examples
#' \dontrun{
#' circ.tree.selectNodes(expr = matrix_clusters, expr_saturated = matrix_clusters01,
#'     cell_clusters = cell_clustering)}

circ.tree.selectNodes <- function(exprs, exprs_saturated, cell_clusters, dist_method = "euclidean", hclust_method = "average"){
  ## Circular hyerarchical clustering tree
  #median expression of each marker for each cell population
  expr_medianL <- data.frame(exprs, cell_clustering = cell_clusters) %>%
    group_by(cell_clustering) %>% summarize_all(list(median)) %>% as.data.frame()
  expr01_medianL <- data.frame(exprs_saturated, cell_clustering = cell_clusters) %>%
    group_by(cell_clustering) %>% summarize_all(list(median)) %>% as.data.frame()
  
  #calculate cluster frequencies (and annotation for heatmap)
  clustering_propL <- data.frame(node = 1:length(levels(cell_clusters)), prop.table(table(cell_clusters))*100)
  
  #hierarchical clustering on clusters
  dL <- dist(expr_medianL, method = dist_method)
  cluster_rowsL <- hclust(dL, method = hclust_method)
  expr_heatL <- as.matrix(expr01_medianL)
  rownames(expr_heatL) <- expr01_medianL$cell_clustering
  
  #hyerarchical tree building
  hca <- as.phylo(cluster_rowsL)
  hca$tip.label <- rownames(expr_heatL)
  
  g <- ggtree(hca, layout = "fan", branch.length = 1) + 
    geom_text2(aes(subset=!isTip, label = node), hjust = -.3) + geom_tiplab()

  ## check X11 active to redirect output
    if("try-error" %in% class(suppressWarnings(try(x11(), silent = T)))){
      cat("X11 is not active, boxplot is saved in -> circTree_selectNodes.", deparse(substitute(cell_clusters)), ".jpg\n", sep = "")
      suppressMessages(ggsave(paste0("circTree_selectNodes.", deparse(substitute(cell_clusters)), ".jpg"), device = "jpeg", plot = g))
    }else{
      print(g)
    }
}


#' circ.tree
#'
#' This function plots a circular dendrogram and a heatmap according expression clusters previously calculated. Every leaf has a different point size regarding the frequency for each cluster.
#' For additional information go to \href{https://guangchuangyu.github.io/software/ggtree/documentation/}{\code{\link{ggtree}}} package.
#' @param exprs A data.frame with marker expression values for each cell
#' @param exprs_saturated A data.frame with marker expression values for each cell where those expression values greater or smaller that 1 and 0, respectively, must be coherced to 1 and 0.
#' @param cell_clusters A vector with clusters identified in a previous clustering method.
#' @param dist_method The distance measure to be used. Possible values are "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param hclust_method The agglomeration method to be used. Possible values are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".
#' @param open_angle Open angle for circular layout. Default = 100
#' @param dendro_labels Logical whether adding labels to dendrogram. Default = \code{FALSE}
#' @param nodes Those nodes limiting different color areas (value can be manually selected from \code{\link[FlowCT:circ.tree]{FlowCT::circ.tree.selectNodes()}}). Default = NULL
#' @param set.seed \code{\link[base:set.seed]{set.seed()}}, for selecting colors for color areas. Default = 333
#' @keywords circular tree
#' @keywords dendrogram
#' @keywords node
#' @export
#' @import dplyr
#' @import stats
#' @import ggtree
#' @importFrom ape as.phylo
#' @examples
#' \dontrun{
#' circ.tree(expr = matrix_clusters, expr_saturated = matrix_clusters01,
#'     cell_clusters = cell_clustering, dendro_labels = TRUE)
#' circ.tree(expr = matrix_clusters, expr_saturated = matrix_clusters01,
#'     cell_clusters = cell_clustering, nodes = c(29,33,35))}

circ.tree <- function(exprs, exprs_saturated, cell_clusters, dist_method = "euclidean", hclust_method = "average", 
                      open_angle = 100, dendro_labels = FALSE, nodes = NULL, set.seed = 333){
  ## Circular hyerarchical clustering tree
  #median expression of each marker for each cell population
  expr_medianL <- data.frame(exprs, cell_clustering = cell_clusters) %>%
    group_by(cell_clustering) %>% summarize_all(list(median)) %>% as.data.frame()
  expr01_medianL <- data.frame(exprs_saturated, cell_clustering = cell_clusters) %>%
    group_by(cell_clustering) %>% summarize_all(list(median)) %>% as.data.frame()
  
  #calculate cluster frequencies (and annotation for heatmap)
  clustering_propL <- data.frame(node = 1:length(levels(cell_clusters)), prop.table(table(cell_clusters))*100)
  
  #hierarchical clustering on clusters
  dL <- dist(expr_medianL, method = dist_method)
  cluster_rowsL <- hclust(dL, method = hclust_method)
  expr_heatL <- as.matrix(expr01_medianL)
  rownames(expr_heatL) <- expr01_medianL$cell_clustering
  
  #hyerarchical tree building
  hca <- as.phylo(cluster_rowsL)
  hca$tip.label <- rownames(expr_heatL)
  
  p1 <- ggtree(hca, layout = "fan", open.angle = open_angle, branch.length = 1) + 
    geom_point(aes(shape = "1", color = cell_clusters, size = Freq)) +
    scale_color_manual(values = colors_palette) + 
    theme(legend.position="right")
  
  p1 <- p1 %<+% clustering_propL #add dataframe for geom_point level
  
  set.seed(set.seed)
  if(!is.null(nodes)){
    for(i in nodes){
      p1 <- p1 + geom_hilight(node = i, fill = sample(colors_palette, 1), alpha = .6)
    }
  }
  
  if(dendro_labels == TRUE){
    p1 + geom_tiplab2(offset=0.1, align = F, size=3)
  }
  
  expr_tree_plot <- expr01_medianL[,-1] #add saturated expression to tree heatmap
  rownames(expr_tree_plot) <- rownames(expr_heatL)
  
  g <- gheatmap(p1, expr_tree_plot, offset = 0.5, width=1, font.size=2, colnames_angle=0, hjust=0,
           colnames_position = "top", high="#b30000", low="#fff7f3") +
    theme(legend.position = NULL)

  ## check X11 active to redirect output
  if("try-error" %in% class(suppressWarnings(try(x11(), silent = T)))){
    cat("X11 is not active, boxplot is saved in -> circTree_heatmap.", deparse(substitute(cell_clusters)), ".jpg\n", sep = "")
    suppressMessages(ggsave(paste0("circTree_heatmap.", deparse(substitute(cell_clusters)), ".jpg"), device = "jpeg", plot = g))
  }else{
    print(g)
  }
}

#' cor.plot.conditions
#'
#' This function draws a correlation plot (based in \code{\link[corrplot:corrplot]{corrplot::corrplot()}}) between specified conditions.
#' For additional information go to \href{https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html}{\code{\link{corrplot}}} package.
#' @param data A data.frame with proportions for each cell type.
#' @param colname_condition Colname of the column containing differential condition.
#' @param colname_patientID Colname of the column containing sample (or patient) identification.
#' @param conditions A vector with those conditions within colname_condition.
#' @param metadata A data.frame with additional information for cell proportions (it must include colname_condition and colname_patientID)
#' @keywords correlation plot
#' @keywords corr
#' @export
#' @importFrom Hmisc rcorr
#' @importFrom corrplot corrplot
#' @import grDevices
#' @examples
#' \dontrun{
#' cor.plot.conditions(data = dataset_final, colname_condition = "condition",
#'     colname_patientID = "patient_id", conditions = c("PB", "BM"), metadata = md)}

cor.plot.conditions <- function(data, colname_condition, colname_patientID, conditions, metadata){
  dataset1 <- data[data[,colname_condition] == conditions[1],]
  dataset2 <- data[data[,colname_condition] == conditions[2],]
  
  dataset_reduced <- merge(dataset1, dataset2, by = colname_patientID, suffixes = paste0(".", conditions))
  no_cols <- paste(c(colnames(metadata), "unclassified"), collapse = "|")
  dataset_reduced <- dataset_reduced[,!grepl(no_cols, colnames(dataset_reduced))]
  
  #build the matrix to perform correlation PB/BM
  corr <- rcorr(scale(dataset_reduced)) #Pearson's correlation
  
  corr_r <- as.matrix(corr$r[grepl(conditions[1], rownames(corr$r)), grepl(conditions[2], colnames(corr$r))]) 
  pval <- as.matrix(corr$P[grepl(conditions[1], rownames(corr$P)), grepl(conditions[2], colnames(corr$P))])
  
  ## check X11 active to redirect output
  if("try-error" %in% class(suppressWarnings(try(x11(), silent = T)))){
    cat("X11 is not active, boxplot is saved in -> corPlot_conditions.", deparse(substitute(cell_clusters)), ".jpg\n", sep = "")
    
    jpeg(paste0("corPlot_conditions.", deparse(substitute(cell_clusters)), ".jpg"), width = 900, heigh = 900, pointsize = 20)
    corrplot(corr_r, order = "hclust", p.mat = pval, insig = "label_sig",
             sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white",
             tl.col="black", tl.cex=.7, tl.offset=0.5, tl.srt=45, 
             col = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                          "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                          "#4393C3", "#2166AC", "#053061")))(100))
    invisible(dev.off())
  }else{
    corrplot(corr_r, order = "hclust", p.mat = pval, insig = "label_sig",
             sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white",
             tl.col="black", tl.cex=.7, tl.offset=0.5, tl.srt=45, 
             col = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                          "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                          "#4393C3", "#2166AC", "#053061")))(100))
  }
}



#' barplot.cell.pops
#'
#' This function calculates cluster proportions for each identified cluster and draws it over a stacked barplot by condition.
#' @param cell_clusters A vector with clusters identified in a previous clustering method.
#' @param metadata A data.frame with additional information for cell proportions (it must include a column with the condition and named as "condition").
#' @param colname_sampleID Colname of the column containing sample identification.
#' @keywords proportions
#' @keywords barplot
#' @export barplot.cell.pops
#' @method barplot cell.pops
#' @import ggplot2
#' @importFrom data.table melt
#' @examples
#' \dontrun{
#' props_final <- barplot.cell.pops(cell_clusters = final_metadata$final_cluster,
#'     metadata = final_metadata[,1:3], colname_sampleID = "sample_id")}

barplot.cell.pops <- function(cell_clusters, metadata, colname_sampleID){
  prop_tableL <- prop.table(table(cell_clusters, metadata[,colname_sampleID]), margin = 2)*100
  
  ggdfL <- data.table::melt(prop_tableL, value.name = "proportion")
  colnames(ggdfL)[2] <- colname_sampleID
  
  mmL <- match(ggdfL[,colname_sampleID], metadata[,colname_sampleID]) #add other infos
  ggdfL <- data.frame(metadata[mmL,], ggdfL)
  
  g <- ggplot(ggdfL, aes_string(x = "sample_id", y = "proportion", fill = "cell_clusters")) +
    geom_bar(stat = "identity") +
    facet_wrap(~ condition, scales = "free_x") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = colors_palette)
  
  ## check X11 active to redirect output
  if("try-error" %in% class(suppressWarnings(try(x11(), silent = T)))){
    cat("X11 is not active, boxplot is saved in -> barPlot_cell_prop.", deparse(substitute(cell_clusters)), ".jpg\n", sep = "")
    suppressMessages(ggsave(paste0("barPlot_cell_prop.", deparse(substitute(cell_clusters)), ".jpg"), device = "jpeg", plot = g))
  }else{
    print(g)
  }
  
  return(prop_tableL)
}


#' cell.count.bx
#'
#' This function draws a boxplot according to the cell count for a specified condition.
#' @param data A data.frame with the metadata for each cell.
#' @param counts_by Condition to perform the counting. Default = "sample_id"
#' @param metadata A data.frame with additional information for cell proportions.
#' @param labels Labels to add to the points on the boxplot, it must be taken from \code{metadata}. Default = \code{NULL}
#' @keywords cell count
#' @keywords boxplot
#' @export
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @examples
#' \dontrun{cell.count.bx(metadata_sc, counts_by = "sample_id", metadata = md)}

cell.count.bx <- function(data, counts_by = "sample_id", metadata, labels = NULL){
  ggdf <- data.frame(metadata, cell_counts = as.numeric(table(data[,counts_by])))
  
  g <- ggplot(ggdf, aes_string(x = "condition", y = "cell_counts", fill = "condition"))+
    geom_boxplot()+
    scale_fill_manual(values = colors_palette, drop = FALSE)+ 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))+
    geom_jitter(shape=16, position=position_jitter(0.2))   
  
  if(!is.null(labels)){
    g + geom_label_repel(aes(label = eval(parse(text = labels))), alpha=0.75, fontface = 'bold', color = 'black')
  }else{
    ## check X11 active to redirect output
    if("try-error" %in% class(suppressWarnings(try(x11(), silent = T)))){
      cat("X11 is not active, boxplot is saved in -> cell_count_boxplot.", deparse(substitute(data)), ".jpg\n", sep = "")
      suppressMessages(ggsave(paste0("cell_count_boxplot.", deparse(substitute(data)), ".jpg"), device = "jpeg", plot = g))
    }else{
      print(g)
    }
  }
}
