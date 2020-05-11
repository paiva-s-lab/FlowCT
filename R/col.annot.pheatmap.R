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