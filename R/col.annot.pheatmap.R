#' Create a palette of colors suitable for doing a heatmap
#' @export

col.annot.pheatmap <- function(metadata, colors = NULL){
  annot_col <- list()
  for(i in colnames(metadata)){
    aux <- as.vector(unlist(unique(metadata[i])))
    if(is.null(colors)) colors2 <- div.colors(length(aux)) else colors2 <- colors
    names(colors2) <- aux
    annot_col[[i]] <- colors2
  }
  return(lapply(annot_col, function(x) x[!is.na(names(x))])) # delete NAs for list with less elements specified in "colors"
}
