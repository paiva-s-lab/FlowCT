#' dr.plotting
#'
#' This function plots the indicated dimensional reduction (DR) from a previously calculated \code{\link[FlowCT:dim.reduction]{FlowCT::dim.reduction()}} object.
#' @param dr A object with DR generated with \code{\link[FlowCT:dim.reduction]{FlowCT::dim.reduction()}} or a \code{data.frame} with DR, expression and metadata information (like the first element list of the object generated with \code{\link[FlowCT:dim.reduction]{FlowCT::dim.reduction()}}).
#' @param dr.calculated String indicating the desired DR to plot (this indicated DR should be prevoulsy calculated to being plotted).
#' @param n.dims Vector indicating the two DR components to plot. Default = \code{c(1,2)}.
#' @param color.by Variable from (from \code{colData(fcs.SE)}) for dots coloring. If \code{color.by = "expression"} (default), plot will be splitted for each marker (\code{facet.by}).
#' @param facet.by Variable from (from \code{colData(fcs.SE)}) for plot spliting. Default = \code{NULL}.
#' @param omit.markers Vector with markers to omit when plotting with \code{color.by = "expression"}. Default = \code{NULL}.
#' @param title Title to add to the plot.
#' @param label.by Variable from (from \code{colData(fcs.SE)}) for dots labeling. Default = \code{NULL}.
#' @param size Point size. Default = \code{0.5}.
#' @param raster Vector indicating if image should be rasterized (logical element) and the number of pixels to consider (numerical element). It is based on (\href{https://github.com/exaexa/scattermore}{scattermore package}). Default = \code{c(T, 1000)}.
#' @keywords dimensional reduction plotting
#' @keywords tSNE
#' @keywords PCA
#' @keywords UMAP
#' @export
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' dr.plotting(dr, dr.calculated = "tSNE", color.by = "condition")
#' dr.plotting(dr, dr.calculated = "UMAP", color.by = "patient_id")
#' dr.plotting(dr, dr.calculated = "PCA", color.by = "SOM", facet.by = "condition")
#' }

dr.plotting <- function(dr, dr.calculated, n.dims = c(1,2), color.by = "expression", facet.by = NULL, omit.markers = NULL, title = "", label.by = NULL, size = 2, raster = c(T, 1000)){
  require(ggplot2)

  if(class(dr) == "list"){
    if(is.null(omit.markers)){
      if(color.by == "expression"){
        dr <- dr$dr_melted
      }else{
        dr <- dr$dr
      }
    }else{
      if(color.by == "expression"){
        dr <- dr$dr_melted[!grepl(paste0(omit.markers, collapse = "|"), dr$dr_melted$antigen),]
      }else{
        dr <- dr$dr[,!grepl(paste0(omit.markers, collapse = "|"), colnames(dr$dr))]
      }
    }
  }else{
    dr <- dr
  }
  
  g <- ggplot(dr, aes_string(x = paste0(dr.calculated, n.dims[1]), y = paste0(dr.calculated, n.dims[2]), 
                             color = color.by)) +
    xlab(paste0(dr.calculated, 1)) + ylab(paste0(dr.calculated, 2)) +
    ggtitle(title) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA))
  
  if(is.factor(dr[,color.by])){
    g <- g + scale_color_manual(values = div.colors(length(unique(dr[,color.by]))), name = color.by) +
      guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
  }else{
    g <- g + scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50), name = color.by)
  }
  
  if(!is.null(facet.by)){
    g <- g + facet_wrap(~ eval(parse(text = facet.by)))
  }else if(color.by == "expression" & is.null(facet.by)){
    g <- g + facet_wrap(~ antigen)
  }

  if(!is.null(label.by)){
    g <- g + geom_text(aes_string(label = label.by), nudge_y = 0.05)
  }

  if(raster[1]){
    g <- g + scattermore::geom_scattermore(pointsize = size, pixels = rep(raster[2], 2)) #devtools::install_github('exaexa/scattermore')
  }else{
    g <- g + geom_point(aes_string(color = color.by), size = size)
  }
  
  return(g)
}


#' circ.tree
#'
#' This function plots a circular dendrogram and a heatmap according a \code{FCS.SE}. Every leaf has a different colored point size regarding the cell type and the frequency for each cluster.
#' For additional information go to \href{https://guangchuangyu.github.io/software/ggtree/documentation/}{\code{\link{ggtree}}} package.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which drawing the circular tree. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param dist.method Distance method measurement to be used. Possible values are "euclidean" (default), "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param hclust.method Hierarchical clustering method to be used. Possible values are "average" (default), "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid".
#' @param nodes If \code{"display"} (default), nodes will be numered. If contains a numeric vector with node numbers, areas defined from these nodes will be differently colored.
#' @param open.angle Angle aperture circular layout. Default = \code{100}.
#' @param dendro.labels Logical whether adding labels to dendrogram. Default = \code{FALSE}.
#' @param scale.size Numerical value indicating how much scale points in the dendogram terminal nodes. Default = \code{10}.
#' @keywords circular tree
#' @keywords dendrogram
#' @keywords nodes
#' @keywords hierachical clustering
#' @export
#' @import dplyr
#' @import phytools
#' @import ggtree
#' @import treeio isTip
#' @import ggplot2
#' @examples
#' \dontrun{
#' # step 1: display all node numbers to select how to coloring areas
#' circ.tree(fcs.SE = fcs_seL, cell.clusters = fcs_seL$SOM_L_named, nodes = "display")
#' 
#' # step 2: color areas indicating node numbers
#' circ.tree(fcs.SE = fcs_seL, cell.clusters = fcs_seL$SOM_L_named, nodes = c(17, 25))
#' }

circ.tree <- function(fcs.SE, assay.i = "normalized", cell.clusters, dist.method = "euclidean", hclust.method = "average", 
                      nodes = "display", open.angle = 100, dendro.labels = FALSE, scale.size = 10){
  require(ggtree)
  require(phytools)

  exprs <- t(assay(fcs.SE, i = assay.i))
  exprs_01 <- exprs.saturate(exprs)
  colors_palette <- div.colors(length(cell.clusters))
  
  ## Circular hyerarchical clustering tree
  #median expression of each marker for each cell population
  expr_medianL <- data.frame(exprs, cell_clustering = cell.clusters) %>%
    group_by(.data$cell_clustering) %>% summarize_all(list(stats::median)) %>% as.data.frame(.data)
  expr01_medianL <- data.frame(exprs_01, cell_clustering = cell.clusters) %>%
    group_by(.data$cell_clustering) %>% summarize_all(list(stats::median)) %>% as.data.frame(.data)
  
  #calculate cluster frequencies (and annotation for heatmap)
  clustering_propL <- data.frame(node = 1:length(levels(cell.clusters)), prop.table(table(cell.clusters))*100)
  
  #hierarchical clustering on clusters
  dL <- stats::dist(expr_medianL, method = dist.method)
  cluster_rowsL <- stats::hclust(dL, method = hclust.method)
  expr_heatL <- as.matrix(expr01_medianL)
  rownames(expr_heatL) <- expr01_medianL$cell_clustering
  
  #hyerarchical tree building
  hca <- ape::as.phylo(cluster_rowsL)
  hca$tip.label <- rownames(expr_heatL)
  
  if(nodes == "display"){
    print(ggtree(hca, layout = "fan", branch.length = 1) + geom_text2(aes_string(subset =! "isTip", label = "node"), hjust = -.3) + geom_tiplab())
  }else{
    p1 <- ggtree(hca, layout = "fan", open.angle = open.angle, branch.length = 1) 
    
    p1 <- p1 %<+% clustering_propL #add dataframe for geom_point level
    
    for(i in nodes){
      p1 <- p1 + geom_hilight(node = i, fill = sample(colors_palette, 1), alpha = .6) +
        geom_point(aes_string(shape = "1", color = "cell_clusters", size = "Freq")) +
        scale_size_area(max_size = scale.size) +
        scale_color_manual(values = colors_palette) + 
        theme(legend.position = "right")
    }
    
    if(dendro.labels){
      p1 + geom_tiplab2(offset=0.1, align = F, size=3)
    }
    
    expr_tree_plot <- expr01_medianL[,-1] #add saturated expression to tree heatmap
    rownames(expr_tree_plot) <- rownames(expr_heatL)
    
    print(gheatmap(p1, expr_tree_plot, offset = 0.5, width = 1, font.size = 2, colnames_angle = 0, hjust = 0,
                   colnames_position = "top", high = "#b30000", low = "#fff7f3") +
            theme(legend.position = NULL))
  }
}


#' cor.plot.conditions
#'
#' It draws a correlation plot (between two specified conditions within a \code{FCS.SE} object.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param condition.column Column name from the \code{colData(fcs.SE)} object which contains condition information.
#' @param return.stats Logical indicating if calculated statistics should be returned in a new variable. Default = \code{FALSE}.
#' @keywords correlation plot
#' @keywords corr
#' @export
#' @importFrom grDevices colorRampPalette
#' @examples
#' \dontrun{
#' corplot.conditions(fcs.SE = fcs_se, cell.clusters = fcs_se$SOM_named, 
#'     condition.column = "condition")
#' }

corplot.conditions <- function(fcs.SE, assay.i = "normalized", cell.clusters, condition.column, return.stats = F){
  prop_table <- as.data.frame.matrix(t(barplot.cell.pops(fcs.SE = fcs.SE, cell.clusters = cell.clusters, count.by = "filename", plot = F)))
  
  prop_table_md <- merge(fcs.SE@metadata$reduced_metada, prop_table, by.x = "filename", by.y = "row.names")
  
  conditions <- as.character(unique(prop_table_md[,condition.column]))
  if(length(conditions) != 2) stop("corplot is onlyprepared for comparing ONLY two conditions")
  
  dataset1 <- prop_table_md[prop_table_md[,condition.column] == conditions[1],colnames(prop_table)]
  colnames(dataset1) <- paste0(colnames(dataset1),  ":", conditions[1])
  rownames(dataset1) <- NULL
  dataset2 <- prop_table_md[prop_table_md[,condition.column] == conditions[2],colnames(prop_table)]
  colnames(dataset2) <- paste0(colnames(dataset2), ":", conditions[2])
  rownames(dataset2) <- NULL
  
  dataset <- merge(dataset1, dataset2, by = "row.names")[-1]
  corr <- Hmisc::rcorr(as.matrix(dataset))
  
  # face conditions
  corr_r <- as.matrix(corr$r[grepl(conditions[1], rownames(corr$r)), grepl(conditions[2], colnames(corr$r))]) 
  pval <- as.matrix(corr$P[grepl(conditions[1], rownames(corr$P)), grepl(conditions[2], colnames(corr$P))])
  
  corrplot::corrplot(corr_r, order = "hclust", p.mat = pval, insig = "label_sig",
                           sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white",
                           tl.col="black", tl.cex=.7, tl.offset=0.5, tl.srt=45, 
                           col = colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                                        "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                                        "#4393C3", "#2166AC", "#053061")))(100))
  if(return.stats) return(corr_r)
}



#' barplot.cell.pops
#'
#' This function calculates cluster proportions (or raw counts) for each identified cluster and plot them on a stacked barplot.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param plot Logical indicating whether plotting stacked barplot. Default = \code{TRUE}.
#' @param count.by Variable name (from \code{colData(fcs.SE)}) for calculating proportions (or counts) and drawing the x-axis in the stacked bar plotting.
#' @param facet.by Variable name (from \code{colData(fcs.SE)}) for splitting the stacked bar plotting. Default = \code{NULL}.
#' @param return.mode String for specifying if final resuls should be proportions ("percentage", default) or raw counts ("counts").
#' @keywords proportions
#' @keywords barplot
#' @export barplot.cell.pops
#' @method barplot cell.pops
#' @import ggplot2
#' @importFrom SummarizedExperiment colData assay
#' @examples
#' \dontrun{
#' prop_table <- barplot.cell.pops(fcs.SE = fcs_se, cell.clusters = fcs_se$SOM_named, 
#'     count.by = "sample_id", facet.by = "condition", 
#'     return.mode = "percentage")
#' counts_table <- barplot.cell.pops(fcs.SE = fcs_se, cell.clusters = fcs_se$SOM_named, 
#'     count.by = "condition", return.mode = "counts")
#' }

barplot.cell.pops <- function(fcs.SE, assay.i = "normalized", cell.clusters, plot = T, count.by, facet.by = NULL, return.mode = "percentage"){
  data <- t(assay(fcs.SE, i = assay.i))
  metadata <- colData(fcs.SE)
  colors_palette <- div.colors(length(unique(cell.clusters)))
  
  counts_table <- table(cell.clusters, metadata[,count.by])
  prop_table <- prop.table(counts_table, margin = 2)*100
  
  if(plot){
    ggdf <- data.table::melt(prop_table, value.name = "proportion")
    colnames(ggdf)[2] <- count.by
    
    mm <- match(ggdf[,count.by], metadata[,count.by]) #add other infos
    ggdf <- data.frame(metadata[mm,], ggdf)
    
    g <- ggplot(ggdf, aes_string(x = count.by, y = "proportion", fill = "cell.clusters")) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(values = colors_palette)
    if(is.null(facet.by)) print(g) else print(g + facet_wrap(~ eval(parse(text = facet.by)), scales = "free_x"))
  }
  if(return.mode == "percentage"){
    return(prop_table)
  }else if(return.mode == "counts"){
    return(counts_table)
  }else{
    cat("Please, specify a valid 'return.mode' value, i.e.: counts or percentage")
  }
}

#' cell.count.bx
#'
#' This function draws a boxplot according to the cell count for a specified condition.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param x.axis Variable name (from \code{colData(fcs.SE)}) for drawing the x-axis in plot.
#' @param color.by Variable name (from \code{colData(fcs.SE)}) for coloring the plot. Default = \code{x.axis} (i.e., the same specfied in \code{x.axis}).
#' @param label.by Variable from (from \code{colData(fcs.SE)}) for labeling. Default = \code{NULL}.
#' @param limits Numeric vector with limits for plotting (minimum, maximum). Default = \code{NULL} (i.e., automatically calculated from data).
#' @keywords cell count
#' @keywords boxplot
#' @export
#' @import ggplot2
#' @importFrom SummarizedExperiment colData assay
#' @examples
#' \dontrun{cell.count.bx(fcs_seL, assay.i = "normalized", x.axis = "condition")}

cell.count.bx <- function(fcs.SE, assay.i = "normalized", x.axis, color.by = x.axis, label.by = NULL, limits = NULL){
  data <- merge(t(assay(fcs.SE, i = assay.i)), colData(fcs.SE), by = "row.names")[,-1]
  ggdf <- data.frame(data[!duplicated(data[,"filename"]),], cell_counts = as.numeric(table(data[,"filename"])))
  
  if(is.null(limits)) limits <- c(min(ggdf$cell_counts), max(ggdf$cell_counts))
  
  g <- ggpubr::ggboxplot(ggdf, x = x.axis, y = "cell_counts", fill = color.by) +
    ggpubr::stat_compare_means(label.x = 1.7, label.y = max(ggdf$cell_counts)) +
    scale_y_continuous(limits = limits) + 
    scale_fill_manual(values = div.colors(length(unique(ggdf[,color.by]))), drop = FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    geom_jitter(shape=16, position=position_jitter(0.2))
  
  if(!is.null(label.by)){
    print(g + ggrepel::geom_label_repel(aes(label = eval(parse(text = label.by))), alpha=0.75, fontface = 'bold', color = 'black'))
  }else{
    print(g)
  }
}

#' multidensity
#'
#' This function draws multidensity plot with all FCS files included in a \code{FCS.SE} object. If there are more files than limit specified in \code{ridgeline.lim}, instead of plottting density in a ridge-way all density lines will be overlapped.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param show.markers Vector with markers to plot. Default = \code{"all"}.
#' @param color.by Variable name (from \code{colData(fcs.SE)}) for lines coloring. 
#' @param subsampling Numeric value indicating how many events use to calculate density lines and speed up plotting. Default = \code{NULL}.
#' @param interactive Logical indicating if the user can interact with the (only overlapping-lines) plot. Default = \code{FALSE}.
#' @param ridgeline.lim Numeric value specifying the limit for shifting from ridgeline-mode to overlapping-lines plot. Default = \code{15}.
#' @keywords marker normalization
#' @keywords marker density
#' @keywords marker alignment
#' @export
#' @import ggplot2
#' @importFrom SummarizedExperiment colData assay
#' @importFrom stats median
#' @examples
#' \dontrun{
#' multidensity(fcs.SE = fcs_se, assay.i = "normalized", subsampling = 1000)
#' multidensity(fcs_se, assay.i = 2, color.by = "file_name", ridgeline.lim = 0, 
#'     show.markers = c("CD62L", "CD4"), interactive = T)
#' }

multidensity <- function(fcs.SE, assay.i, show.markers = "all", color.by = NULL, subsampling = NULL, interactive = F, ridgeline.lim = 15){
  if(show.markers == "all") show.markers <- rownames(fcs.SE)
  if(!is.null(subsampling)) suppressMessages(fcs.SE <- sub.samples(fcs.SE, subsampling = subsampling))
  
  data <- t(assay(fcs.SE, i = assay.i))
  data2 <- merge(data, colData(fcs.SE), by = "row.names")[,-1]
  
  # prepare tables: for plotting and with median values for each marker
  median_df <- data.frame(antigen = show.markers, median = apply(data[,show.markers], 2, median))
  ggdf <- data.table::melt(data2, measure.vars = show.markers, value.name = "expression", variable.name = "antigen")
  
  if(length(unique(fcs.SE$filename)) > ridgeline.lim){
    g <- ggplot(data = ggdf[grepl(paste0(show.markers, collapse = "|"), ggdf$antigen),], 
                aes_string(x = "expression", color = color.by, group = "filename")) + 
      # geom_density(size = 0.5) +
      stat_density(geom = "line", position = "identity", size = 0.5) +
      facet_wrap(~ antigen, scales = "free") +
      geom_vline(data = median_df, aes(xintercept = median), linetype = 2, color = "gray55") +
      scale_color_manual(name = color.by, values = div.colors(unique(length(ggdf[,color.by])))) +
      theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                              strip.text = element_text(size = 7), axis.text = element_text(size = 5))
    
    if(interactive) plotly::ggplotly(g) else print(g)
  }else{
    suppressMessages(print(ggplot(ggdf, aes_string(x = "expression", y = "filename")) + 
                             ggridges::geom_density_ridges(alpha = 0.7) +
                             facet_wrap(~ antigen, scales = "free") +
                             geom_vline(data = median_df, aes(xintercept = median), linetype = 2, color = "gray55") +
                             theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                                                     strip.text = element_text(size = 7), axis.text = element_text(size = 7))))
  }
}


#' median.heatmap
#'
#' This function draws a heatmap with median values for each FCS file or for identified cluster with \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}}
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed). Default = \code{NULL} (i.e., median values will be calculated for each FCS file).
#' @param markers.to.use Vector with markers to use. Default = \code{"all"}.
#' @param not.metadata Vector with variable names (from \code{colData(fcs.SE)}) for not including in the heatmap annotation. Default = \code{"filename"}.
#' @keywords heatmap
#' @keywords cell cluster percentages
#' @keywords median expression values
#' @export median.heatmap
#' @import dplyr
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' \dontrun{
#' # option 1: general heatmap (by FCS file)
#' median.heatmap(fcs_se, not.metadata = c("sample_id", "file_name"))
#' 
#' # option 2: heatmap by SOM-identified clusters
#' median.heatmap(fcs.SE = fcs_se, assay.i = "normalized", cell.clusters = fcs_se$SOM)
#' }

median.heatmap <- function(fcs.SE, assay.i = "normalized", cell.clusters = NULL, markers.to.use = "all", not.metadata = "filename"){
  data <- t(assay(fcs.SE, i = assay.i))
  metadata <- fcs.SE@metadata$reduced_metadata
  
  if(markers.to.use == "all") markers.to.use <- colnames(data)
  
  ## prepare median tables
  if(is.null(cell.clusters)){
    med <- median.values(fcs.SE, assay.i = assay.i)
  }else{
    expr_median <- data.frame(cell_clusters = cell.clusters, data[,markers.to.use]) %>%
      group_by(.data$cell_clusters) %>% summarize_all(list(stats::median)) %>% as.data.frame(.data)
    
    expr_saturated_median <- data.frame(cell_clusters = cell.clusters, exprs.saturate(data[,markers.to.use])) %>%
      group_by(.data$cell_clusters) %>% summarize_all(list(stats::median)) %>% as.data.frame(.data)
  }
  
  ## heatmap
  if(is.null(cell.clusters)){
    annotation_colors <- col.annot.pheatmap(metadata[,!(colnames(metadata) %in% not.metadata), drop = F])
    color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
    
    print(pheatmap::pheatmap(t(med[,markers.to.use]), color = color, display_numbers = FALSE,
                   number_color = "black", fontsize_number = 5, clustering_method = "average",
                   annotation = metadata[,!(colnames(metadata) %in% not.metadata), drop = F], 
                   annotation_colors = annotation_colors, 
                   show_colnames = F))
  }else{
    ## calculate cluster frequencies
    clustering_table <- table(cell.clusters)
    clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
    labels_row <- paste0(expr_saturated_median$cell_clusters, " (", clustering_prop ,"%)")
    
    d <- stats::dist(expr_median[,markers.to.use], method = "euclidean")
    cluster_rows <- stats::hclust(d, method = "average")
    expr_heat <- as.matrix(expr_saturated_median[,markers.to.use])
    # rownames(expr_heat) <- paste0("c.", expr_saturated_median$cell_clusters) #force rownames to not crash heatmap (??)
    rownames(expr_heat) <- rownames(expr_saturated_median) #force rownames to not crash heatmap (??)
    
    ## annot colors
    annot_row <- expr_saturated_median[,"cell_clusters", drop = F]
    # rownames(annot_row) <- paste0("c.", rownames(annot_row)) #force rownames to not crash heatmap (??)
    annotation_colors <- col.annot.pheatmap(expr_saturated_median[,"cell_clusters", drop = F])
    color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
    legend_breaks <- seq(from = 0, to = 1, by = 0.1)
    
    print(pheatmap::pheatmap(expr_heat, color = color, annotation_legend = F,
                   cluster_cols = FALSE, cluster_rows = cluster_rows, labels_row = labels_row,
                   display_numbers = FALSE, number_color = "black",
                   fontsize = 8, fontsize_number = 6, legend_breaks = legend_breaks,
                   annotation_row = annot_row, annotation_colors = annotation_colors))
  }
}


# 'median.dr
#'
#' This function draws a plot from internally calculated dimensional reduction (DR).
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param markers.to.use Vector with markers to use. Default = \code{"all"}.
#' @param dr.method DR method to calculate. Possible values are "PCA", "tSNE" or "UMAP", individually, or "all". Take into account the number of samples for this DR, for few samples t-SNE and UMAP are useless.
#' @param perplexity.tsne Value for perplexity parameter in tSNE calculation (\href{https://distill.pub/2016/misread-tsne/}{more information}). For little number of samples this DR (t-SNE) cannot be calculated because preplexity restrictions. Default = \code{15}.
#' @param n.neighbors.umap Value for neighbors parameter in UMAP calculation (\href{https://www.math.upenn.edu/~jhansen/2018/05/04/UMAP/}{more information}). Default = \code{50}.
#' @param color.by Variable name (from \code{colData(fcs.SE)}) for lines coloring. 
#' @param label.by Variable name (from \code{colData(fcs.SE)}) for labelling. 
#' @param size Point size. Default = \code{3}.
#' @keywords dimensional reduction
#' @keywords median expression values
#' @export median.dr
#' @examples
#' \dontrun{
#' median.dr(fcs_se, color.by = "condition")
#' }

median.dr <- function(fcs.SE, assay.i = "normalized", markers.to.use = "all", dr.method = "PCA", perplexity.tsne = 15, n.neighbors.umap = 50, color.by, label.by = NULL, size = 3){
  data <- t(SummarizedExperiment::assay(fcs.SE, i = assay.i))
  metadata <- fcs.SE@metadata$reduced_metadata
  
  if(markers.to.use == "all") markers.to.use <- colnames(data)
  
  ## prepare median tables
  med <- median.values(fcs.SE, assay.i = assay.i)
  
  ## dr
  if(!is.null(dr.method)){
    if(grepl("PCA|tSNE|UMAP|all", dr.method)){
      dr <- dim.reduction(med, metadata = metadata, 
                          dr.method = dr.method, markers.to.use = markers.to.use, 
                          perplexity.tsne = perplexity.tsne, 
                          n.neighbors.umap = n.neighbors.umap)
      print(dr.plotting(dr, dr.calculated = dr.method, color.by = color.by, label.by = label.by, 
                        size = size, raster = F))
    }else{
      stop("Please, specify a valid reduction method: PCA, tSNE or UMAP.")
    }
  }
}


#' col.annot.pheatmap
#' @export

col.annot.pheatmap <- function(metadata, colors.palette = NULL){
  annot_col <- list()
  for(i in colnames(metadata)){
    aux <- as.vector(unlist(unique(metadata[i])))
    if(is.null(colors.palette)){
      colors.palette2 <- c()
      colors.palette2 <- div.colors(length(aux))
    }
    names(colors.palette2) <- aux
    annot_col[[colnames(metadata[i])]] <- colors.palette2
  }
  return(annot_col)
}


# 'sc.Heatmap
#'
#' This function draws a heatmap with single-cell fluorescence values.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param markers.to.use Vector with markers to use. Default = \code{"all"}.
#' @param not.metadata Vector with variable names (from \code{colData(fcs.SE)}) for not including in the heatmap annotation. Default = \code{c("filename", "cell_id", "sample_id")}.
#' @param clustering.method Clustering method for rows and columns clustering within the heatmap. Possible values are "average" (default), "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid".
#' @param subsampling Numeric value indicating how many events use to draw heatmap and speed up plotting. Default = \code{100}.
#' @param color Color vector for coloring expression values within the heatmap. Default = \code{NULL} (i.e., scale "YlGnBu" from \code{RColorBrewer}).
#' @keywords single-cell expression
#' @keywords heatmap
#' @export
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom SummarizedExperiment colData assay
#' @examples
#' \dontrun{
#' sc.heatmap(fcs_se100)
#' }

sc.heatmap <- function(fcs.SE, assay.i = "normalized", markers.to.use = "all", not.metadata = c("filename", "cell_id", "sample_id"), 
                       clustering.method = "average", subsampling = 100, color = NULL){
  if(!is.null(subsampling)) suppressMessages(fcs.SE <- sub.samples(fcs.SE, subsampling = subsampling))
  if(markers.to.use == "all") markers.to.use <- rownames(fcs.SE)
  if(is.null(color)) color <- colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100)
  data <- t(assay(fcs.SE, i = assay.i))
  metadata <- as.data.frame(colData(fcs.SE))
    
  annot_col <- col.annot.pheatmap(metadata[,!(colnames(metadata) %in% not.metadata)])
  
  pheatmap::pheatmap(t(data[,markers.to.use]), 
           annotation_col = metadata[,!(colnames(metadata) %in% not.metadata)], 
           annotation_colors = annot_col, color = color, display_numbers = FALSE, 
           number_color = "black", fontsize_number = 5, 
           show_colnames = F, clustering_method = clustering.method, 
           treeheight_row = 0, treeheight_col = 0)

}


# 'boxplot.cell.clustering
#'
#' It draws a boxplot with cell clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}}.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param condition.column Column name from the \code{colData(fcs.SE)} object which contains condition information.
#' @param pvalue.cutoffs List of P-value cutoffs and their symbols for indicante significances within the plot. Default = \code{list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("xxxx", "***", "**", "*", "ns"))}.
#' @param color.by Variable name (from \code{colData(fcs.SE)}) for lines coloring. 
#' @param geom.point Logical indicating if adding points to boxplot. Default = \code{TRUE}.
#' @param facet Logical indicating if splitting boxplots by cell clusters. Default = \code{FALSE}.
#' @param facet.free.scale If \code{facet = TRUE}, string indicating how scales would be shared across all facets. Possible values: "free_x" (default), "free_y" and "free".
#' @param shape.by Variable name (from \code{colData(fcs.SE)}) for dot shaping. Default = \code{NULL}.
#' @param y.limits Numeric vector with limits for y-axis (minimum, maximum). Default = \code{NULL}.
#' @param show.stats Significances should be added to boxplots? Default = \code{TRUE}.
#' @param return.stats Logical indicating if calculated statistics should be returned in a new variable. Default = \code{TRUE}.
#' @param plot.only.sig Vector indicating if only significant cell clusters should be displayed (logical element) and the P-value cutoff for selecting those ones (numerical element). Default = \code{c(F, 0.05)}.
#' @keywords differential boxplot
#' @keywords cell clusters distributions
#' @export boxplot.cell.clustering
#' @import ggplot2
#' @examples
#' \dontrun{
#' # option 1: show all cell clusters and return statistics
#' bx_sig <- boxplot.cell.clustering(fcs.SE = fcs_se, cell.clusters = fcs_se$SOM_named, facet = T, 
#'     return.stats = T, 
#'     facet.free.scale = "free")
#' 
#' # option 2: show only those significant cell clusters
#' boxplot.cell.clustering(fcs.SE = fcs_se, cell.clusters = fcs_se$SOM_named, return.stats = F, 
#'     plot.only.sig = c(T, 0.1))
#' }

boxplot.cell.clustering <- function(fcs.SE, assay.i = "normalized", cell.clusters, condition.column = "condition", 
                                    pvalue.cutoffs = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("xxxx", "***", "**", "*", "ns")), 
                                    color.by = condition.column, geom.point = T, facet = F, facet.free.scale = "free_x", shape.by = NULL, y.limits = NULL,
                                    show.stats = T, return.stats = T, plot.only.sig = c(F, 0.05)){
  ## prepare tables
  prop_table <- as.data.frame.matrix(t(barplot.cell.pops(fcs.SE = fcs.SE, cell.clusters = cell.clusters, count.by = "filename", plot = F)))
  
  prop_table_md <- merge(fcs.SE@metadata$reduced_metada, prop_table, by.x = "filename", by.y = "row.names")
  
  prop_table_mdm <- data.table::melt(prop_table_md, id.vars = colnames(fcs.SE@metadata$reduced_metada), 
                         variable.name = "cell_cluster", value.name = "proportion")
  
  ## statistics table
  resultskw <- matrixTests::col_kruskalwallis(prop_table_md[,as.character(unique(cell.clusters))], prop_table_md[,condition.column])
  KWsig <- rownames(resultskw[resultskw$pvalue < plot.only.sig[2],])
  
  if(length(unique(prop_table_md[,condition.column])) > 2){
    kw_posthoc <- list()
    for(i in KWsig) kw_posthoc[[i]] <- stats::pairwise.wilcox.test(prop_table_md[,i], prop_table_md[,condition.column])
  }
  
  ## plotting
  colors_palette <- div.colors(length(unique(prop_table_md[,condition.column])), set.seed = 3)
  
  if(plot.only.sig[1]){
    prop_table_mdm <- prop_table_mdm[prop_table_mdm$cell_cluster %in% KWsig,]
  }
  
  g <- ggplot(prop_table_mdm, aes_string("cell_cluster", "proportion", color = color.by, fill = color.by)) + 
    geom_boxplot(alpha = 0.6) + labs(x = "cell clusters", y = "Proportion") +
    theme_bw()
  
  if(!is.null(y.limits)) g <- g + scale_y_continuous(limits = c(y.limits))
  
  if(geom.point) g <- g + geom_point(aes_string("cell_cluster", "proportion", color = color.by), 
                                     alpha = 0.8, position = position_jitterdodge())
  
  if(!is.null(shape.by)) g <- g + geom_point(aes_string("cell_cluster", "proportion", color = color.by,
                                                        shape = shape.by), alpha = 0.8, position = position_jitterdodge())
  
  if(facet){
    g <- g + facet_wrap(~ cell_cluster, scales = facet.free.scale) +
                  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  }else{
    g <- g + coord_flip()
  } 
  
  if(show.stats) g <- g + ggpubr::stat_compare_means(label = "p.signif", symnum.args = pvalue.cutoffs)
  
  print(g + scale_fill_manual(values = colors_palette) + scale_color_manual(values = colors_palette))
  
  if(return.stats & length(unique(prop_table_md[,condition.column])) > 2){
    return(list(kruskal_results = resultskw, kw_posthoc_results = kw_posthoc))
  }else if(return.stats){
    return(resultskw)
  }
}

# 'dumbPlot.cell.clustering
#'
#' It draws a Dumbbell plot according condition for each cell cluster identified.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param condition.column Column name from the \code{colData(fcs.SE)} object which contains condition information. De
#' @param psig.cutoff P-value cutoff. Default = \code{0.05}.
#' @param return.stats Logical indicating if calculated statistics should be returned in a new variable. Default = \code{TRUE}.
#' @keywords differential dotplot
#' @keywords Dumbbell plot
#' @keywords longitudinal dotplot
#' @export
#' @import ggplot2
#' @examples
#' \dontrun{
#' diffDots.cell.clustering(fcs.SE = fcs_se, cell.clusters = fcs_se$SOM_named, return.stats = F)
#' }

dumbPlot.cell.clustering <- function(fcs.SE, assay.i = "normalized", cell.clusters, condition.column, psig.cutoff = 0.05, return.stats = T){
  ## prepare tables
  prop_table <- as.data.frame.matrix(t(barplot.cell.pops(fcs.SE, cell.clusters, count.by = "filename", plot = F, assay.i = assay.i)))
  
  prop_table_md <- merge(fcs.SE@metadata$reduced_metada, prop_table, by.x = "filename", by.y = "row.names")
  
  dfm <- data.table::melt(prop_table_md, measure.vars = unique(cell.clusters))
  dfma <- stats::aggregate(dfm$value ~ dfm$variable + dfm[,condition.column], FUN = stats::median)
  colnames(dfma) <- c("variable", "condition", "value")
  dfma <- transform(dfma, pct = log(stats::ave(dfma$value, dfma[,condition.column], FUN = function(x) x/sum(x)*100))) #transform to percentaje
  
  ## statistics table
  resultskw <- matrixTests::col_kruskalwallis(prop_table_md[,as.character(unique(cell.clusters))], prop_table_md[,condition.column])
  KWsig <- rownames(resultskw[resultskw$pvalue < psig.cutoff,])
  dfma$sig <- ifelse(dfma$variable %in% KWsig, "1", "0")
  
  if(length(unique(prop_table_md[,condition.column])) > 2){
    kw_posthoc <- list()
    for(i in KWsig) kw_posthoc[[i]] <- stats::pairwise.wilcox.test(prop_table_md[,i], prop_table_md[,condition.column])
  }
  
  # beta: calculating post-hoc for multiple conditions and draw colored lines between each condition and not general
  # KW_ph_sig <- c()
  # if(length(unique(prop_table_md[,"condition"])) > 2){
  #   for(i in names(kw_posthoc)){
  #     aux_stats <- melt(kw_posthoc[[i]]$p.value)
  #     aux_stats <- aux_stats[aux_stats$value < 0.1,]
  #   }
  # }
  
  ## plotting
  print(ggplot(dfma, aes_string(x = "pct", y = "variable", fill = condition.column, color = "sig")) + 
          geom_line(aes_string(group = "variable"), size = 1) +
          scale_color_manual(values = c("gray63", "brown1"), labels = c("no sig.", "sig.")) +
          geom_point(size = 4, shape = 21, color = "white") +
          scale_fill_manual(values = div.colors(length(unique(prop_table_md[,condition.column])), set.seed = 3), 
                            labels = unique(prop_table_md[,condition.column])) +
          # guides(color = guide_legend(ncol = 1)) + #display legend in one-column format
          theme(panel.background = element_blank(), legend.key = element_blank(),
                panel.grid.major.x = element_line(colour = "gray73", linetype = "dashed", size = 0.3),
                panel.border = element_rect(color = "gray73", fill = NA), 
                axis.text = element_text(color = "gray20", face = "bold"), axis.title = element_text(face = "bold")) +
          xlab("\n Population % (log-transformed)") + ylab("Cell clusters\n") + labs(color = ""))
  
  if(return.stats & length(unique(prop_table_md[,condition.column])) > 2){
    return(list(kruskal_results = resultskw, kw_posthoc_results = kw_posthoc))
  }else if(return.stats){
    return(resultskw)
  }
}


# 'dotplot.DE
#'
#' It draws dotplot for each cell cluster identified and each marker to facilitate identification of cell populations in each cell cluster.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param clusters.named Column name from the \code{initial.fcs.SE} object which contains renamed clusters (through \code{\link[FlowCT:clusters.rename]{FlowCT::clusters.rename()}}).
#' @param markers.to.use Vector with markers to use. Default = \code{"all"}.
#' @param psig.cutoff P-value cutoff. Default = \code{0.05}.
#' @param return.stats Logical indicating if calculated statistics should be returned in a new variable. Default = \code{FALSE}.
#' @param scale.size Numerical value indicating how much scale points. Default = \code{9}.
#' @keywords differential dotplot
#' @keywords cell cluster identification
#' @keywords cell cluster markers
#' @export
#' @import ggplot2
#' @examples
#' \dontrun{
#' de.dotplot(fcs.SE = fcs_se, markers.to.use = c("CD8", "CD27", "CCR4", "CD45RA", "CD4"), 
#'     clusters.named = "SOM_named")
#' }

dotplot.DE <- function(fcs.SE, assay.i = "normalized", clusters.named = "SOM_named", markers.to.use = "all", psig.cutoff = 0.05, return.stats = F, scale.size = 9){
  if(markers.to.use == "all") markers.to.use <- rownames(fcs.SE)
  
  dt <- data.frame()
  for(i in unique(fcs.SE[[clusters.named]])){
    aux_se <- fcs.SE[,fcs.SE[[clusters.named]] == i]
    dt <- rbind(dt, data.frame(t(assay(aux_se, i = assay.i)[markers.to.use,]), pop = i))
  }
  dtm <- data.table::melt(dt)
  
  dtm2 <- data.frame()
  for(i in unique(dtm$variable)){
    s1 <- as.data.frame(summary(stats::lm(value ~ 0 + pop, dtm[dtm$variable == i,]))$coefficients)
    s1$pop <- gsub("pop", "", rownames(s1))
    s1$`t value` <- ifelse(s1$`Pr(>|t|)` < psig.cutoff, s1$`t value`, NA) #not consider those non-significant
    s1$mixstats <- s1$Estimate*s1$`t value`

    
    s2 <- stats::aggregate(value ~ ., dtm[dtm$variable == i,], FUN = stats::median) #collapse to median cell values
    s2 <- merge(s2, s1, by = "pop")
    dtm2 <- rbind(dtm2, s2)
    dtm2$pop <- factor(dtm2$pop)
  }
  
  print(ggplot(dtm2, aes_string("variable", y = "pop", size = "value", color = "mixstats")) + 
    geom_point() + 
    scale_color_continuous(na.value = "gray70", name = "Marker\nimportance") +
    scale_size_area(max_size = scale.size, name = "Median\nfluorescence") +
    labs(x = "", y = "") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  
  if(return.stats) return(dtm2)
}


# 'diffdots.cell.clustering
#'
#' It draws a differential dot plot (longitudinally) according condition for each cell cluster identified.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param condition.column Column name from the \code{colData(fcs.SE)} object which contains condition information. De
#' @param psig.cutoff P-value cutoff. Default = \code{0.05}.
#' @param return.stats Logical indicating if calculated statistics should be returned in a new variable. Default = \code{TRUE}.
#' @keywords differential dotplot
#' @keywords Dumbbell plot
#' @keywords longitudinal dotplot
#' @export
#' @import ggplot2
#' @examples
#' \dontrun{
#' diffDots.cell.clustering(fcs.SE = fcs_se, cell.clusters = fcs_se$SOM_named, return.stats = F)
#' }

diffdots.cell.clustering <- function(fcs.SE, assay.i = "normalized", cell.clusters, condition.column, psig.cutoff = 0.05, return.stats = T){
  ## prepare tables
  prop_table <- as.data.frame.matrix(t(barplot.cell.pops(fcs.SE, cell.clusters, count.by = "filename", plot = F, assay.i = "normalized")))
  
  prop_table_md <- merge(fcs.SE@metadata$reduced_metada, prop_table, by.x = "filename", by.y = "row.names")
  
  dfm <- data.table::melt(prop_table_md, measure.vars = unique(cell.clusters))
  dfma <- stats::aggregate(dfm$value ~ dfm$variable + dfm[,condition.column], FUN = stats::median)
  colnames(dfma) <- c("variable", condition.column, "value")
  dfma <- transform(dfma, pct = log(stats::ave(dfma$value, dfma[,condition.column], FUN = function(x) x/sum(x)*100))) #transform to percentaje
  
  conditions <- unique(dfma[,condition.column])
  
  ## statistics table
  resultskw <- matrixTests::col_kruskalwallis(prop_table_md[,as.character(unique(cell.clusters))], prop_table_md[,condition.column])
  KWsig <- rownames(resultskw[resultskw$pvalue < psig.cutoff,])
  dfma$sig <- ifelse(dfma$variable %in% KWsig, "1", "0")
  
  if(length(unique(prop_table_md[,condition.column])) > 2){
    kw_posthoc <- list()
    for(i in KWsig) kw_posthoc[[i]] <- stats::pairwise.wilcox.test(prop_table_md[,i], prop_table_md[,condition.column])
  }
  
  
  ## plotting
  print(ggplot(dfma, aes_string(x = "condition", y = "pct", fill = "variable")) + 
          geom_point(size = 3, shape = 21, color = "gray63") + 
          # scale_color_manual(values = div.colors(30)) + 
          geom_line(aes_string(group = "variable", color = "sig"), size = 1) + 
          scale_colour_manual(values = c("gray63", "brown1"), labels = c("no sig.", "sig.")) + 
          # geom_label(data = subset(dfma, dfma$sig == 1 & dfma$condition == conditions[1]), 
          #            aes(x = 1, y = pct, label = variable, fill = variable), nudge_x = -0.1) +
          ggrepel::geom_label_repel(data = subset(dfma, dfma$sig == 1 & dfma$condition == conditions[1]),
                     aes_string(x = 1, y = "pct", label = "variable", fill = "variable"), 
                     nudge_x = -0.1, show.legend = F) +
          theme(panel.background = element_blank(), axis.line = element_line(color = "black")) + 
          labs(x = "\nCondition", y = "% of cells (log-transf.)\n", color = "", fill = "Cell clusters"))
  
  if(return.stats & length(unique(prop_table_md[,condition.column])) > 2){
    return(list(kruskal_results = resultskw, kw_posthoc_results = kw_posthoc))
  }else if(return.stats){
    return(resultskw)
  }
}
