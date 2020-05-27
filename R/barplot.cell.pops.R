#' barplot.cell.pops
#'
#' This function calculates cluster proportions (or raw counts) for each identified cluster and plot them on a stacked barplot.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters A vector with clusters identified through \code{\link[FlowCT.v2:fsom.clustering]{FlowCT.v2::fsom.clustering()}} (and, normaly, later renamed).
#' @param plot Logical indicating whether plotting stacked barplot. Default = \code{TRUE}.
#' @param count.by Variable name (from \code{colData(fcs.SCE)}) for calculating proportions (or counts) and drawing the x-axis in the stacked bar plotting.
#' @param facet.by Variable name (from \code{colData(fcs.SCE)}) for splitting the stacked bar plotting. Default = \code{NULL}.
#' @param return.mode String for specifying if final resuls should be proportions ("percentage", default) or raw counts ("counts").
#' @param colors Vector with colors for plotting. Default = \code{NULL} (i.e., it will choose automatically a vector of colors according to \code{\link[FlowCT.v2:div.colors]{FlowCT.v2::div.colors()}}).
#' @keywords proportions
#' @keywords barplot
#' @export barplot.cell.pops
#' @method barplot cell.pops
#' @import ggplot2
#' @importFrom SummarizedExperiment colData assay
#' @importFrom data.table melt as.data.table
#' @examples
#' \dontrun{
#' prop_table <- barplot.cell.pops(fcs.SCE = fcs, cell.clusters = fcs$SOM_named, 
#'     count.by = "sample_id", facet.by = "condition", 
#'     return.mode = "percentage")
#' counts_table <- barplot.cell.pops(fcs.SCE = fcs, cell.clusters = fcs$SOM_named, 
#'     count.by = "condition", return.mode = "counts")
#' }

barplot.cell.pops <- function(fcs.SCE, assay.i = "normalized", cell.clusters, plot = T, count.by, facet.by = NULL, return.mode = "percentage", colors = NULL){
  data <- t(assay(fcs.SCE, i = assay.i))
  metadata <- colData(fcs.SCE)
  if(is.null(colors)) colors <- div.colors(length(unique(cell.clusters)))
  
  counts_table <- table(cell.clusters, metadata[,count.by])
  ggdf <- as.data.frame(as.data.table(prop.table(counts_table, margin = 2)*100))
  colnames(ggdf) <- c("cell.clusters", count.by, "proportion")
  
  if(plot){
    mm <- match(ggdf[,count.by], metadata[,count.by]) #add other infos
    ggdf <- data.frame(metadata[mm,], ggdf)
    
    g <- ggplot(ggdf, aes_string(x = count.by, y = "proportion", fill = "cell.clusters")) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(values = colors)
    if(is.null(facet.by)) print(g) else print(g + facet_wrap(~ eval(parse(text = facet.by)), scales = "free_x"))
  }
  
  if(return.mode == "percentage"){
    return(prop.table(counts_table, margin = 2)*100)
  }else if(return.mode == "counts"){
    return(counts_table)
  }else{
    cat("Please, specify a valid 'return.mode' value, i.e.: counts or percentage")
  }
}
