#' Stacked barplot with cell clusters
#'
#' This function calculates cluster proportions (or raw counts) for each identified cluster and plot them on a stacked barplot.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} (and, normaly, later renamed).
#' @param plot Logical indicating whether plotting stacked barplot. Default = \code{TRUE}.
#' @param count.by Variable name (from \code{colData(fcs.SCE)}) for calculating proportions (or counts) and drawing the x-axis in the stacked bar plotting.
#' @param facet.by Variable name (from \code{colData(fcs.SCE)}) for splitting the stacked bar plotting. Default = \code{NULL}.
#' @param return.mode String for specifying if final resuls should be proportions ("percentage") or raw counts ("counts"). Default = \code{NULL}.
#' @param colors Vector with colors for plotting. Default = \code{NULL} (i.e., it will choose automatically a vector of colors according to \code{\link[FlowCT:div.colors]{FlowCT::div.colors()}}).
#' @keywords proportions
#' @keywords barplot
#' @export barplot.cell.pops
#' @method barplot cell.pops
#' @import ggplot2
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay
#' @examples
#' \dontrun{
#' prop_table <- barplot.cell.pops(fcs.SCE = fcs, cell.clusters = fcs$SOM_named, 
#'     count.by = "sample_id", facet.by = "condition", 
#'     return.mode = "percentage")
#' counts_table <- barplot.cell.pops(fcs.SCE = fcs, cell.clusters = fcs$SOM_named, 
#'     count.by = "condition", return.mode = "counts")
#' }

barplot.cell.pops <- function(fcs.SCE, assay.i = "normalized", cell.clusters, plot = T, count.by, facet.by = NULL, return.mode, colors = NULL){
  data <- t(assay(fcs.SCE, i = assay.i))
  metadata <- colData(fcs.SCE)
  
  counts_table <- table(metadata[,cell.clusters], metadata[,count.by])

  if(plot){
  if(is.null(colors)) colors <- div.colors(length(unique(metadata[,cell.clusters])))

  ggdf <- as.data.frame(metadata) %>% 
  select(count.by, facet.by, cell.clusters) %>% 
  add_count(x = ., eval(parse(text = cell.clusters)), eval(parse(text = count.by)))

  ggdf <- aggregate(n ~ ., ggdf[,-c(ncol(ggdf)-1, ncol(ggdf)-2)], FUN = unique)

    g <- ggplot(ggdf, aes_string(x = count.by, y = "n", fill = cell.clusters)) +
      geom_bar(stat = "identity", position = "fill") +
      theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(values = colors)

    if(!is.null(facet.by)) g <- g + facet_wrap(~ eval(parse(text = facet.by)), scales = "free_x")

    print(g)
  }
  
  if(!missing(return.mode)){
    if(return.mode == "percentage"){
      return(prop.table(counts_table, margin = 2)*100)
    }else if(return.mode == "counts"){
      return(counts_table)
    }else{
      cat("Please, indicate a valid return mode: 'percentage' or 'counts'")
    } 
  }
}
