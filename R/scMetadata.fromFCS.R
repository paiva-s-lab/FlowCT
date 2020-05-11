#' scMetadata.fromFCS
#'
#' It creates, from an initial metadata provided by the user, a global metadata for each event in a (\href{https://www.rdocumentation.org/packages/flowCore/versions/1.38.2/topics/flowSet-class}{flowset}) object.
#' @param flowset A flowset object generated following the (\ref{https://bioconductor.org/packages/devel/bioc/vignettes/flowCore/inst/doc/HowTo-flowCore.pdf}\code{flowCore}{tutorial}).
#' @param metadata Metadata information for each read FCS file.
#' @param add.exprs Logical indicating if expression data must be added to each cell's metadata. Default = \code{TRUE}.
#' @keywords single-cell metadata
#' @keywords flowset metadata
#' @export
#' @examples
#' \dontrun{
#'  sc_metadata <- scMetadata.fromFCS(fcs, metadata, add.exprs = T)
#' }

scMetadata.fromFCS <- function(flowset, metadata, add.exprs = T){
  expr <- flowCore::fsApply(flowset, flowCore::exprs)
  mtd <- data.frame()
  pb <- progress::progress_bar$new(total = length(flowset), format = "Generating single cell metadata [:bar] :percent")
  
  file_count <- 1
  for(i in flowCore::sampleNames(flowset)){
    pb$tick()
    
    aux_md <- subset(metadata, metadata$filename == i)
    mtd_aux <- aux_md[rep(1,nrow(flowset[[i]])),]
    mtd_aux$cell_id <- paste0(file_count, ".", seq_len(nrow(mtd_aux)))
    rownames(mtd_aux) <- mtd_aux$cell_id
    
    mtd <- rbind(mtd, mtd_aux)
    file_count <- file_count + 1
    Sys.sleep(1/10)
  }
  if(add.exprs) mtd_exprs <- cbind(mtd, expr) else mtd_exprs <- mtd
  return(mtd_exprs)
}
