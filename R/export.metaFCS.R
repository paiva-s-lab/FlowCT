#' export.metaFCS
#'
#' It creates a FCS file containing all analysis incorpored to colData(fcs.SE) as well as the dimensional reduction coordinates calculated.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param dr.object Object created with \code{\link[FlowCT:dim.reduction]{FlowCT::dim.reduction()}} function or a table combining DR and experimental metadata information.
#' @param output.name Name for generated FCS file. Important, add the final extension ".fcs".
#' @keywords FCS generation
#' @keywords final FCS
#' @keywords exporting
#' @export
#' @importFrom SummarizedExperiment colData assay
#' @examples
#' \dontrun{
#' export.metaFCS(fcs.SE = fcs_se, dr.object = dr, output.name = "final_file.fcs")
#' }

export.metaFCS <- function(fcs.SE, dr.object, output.name){
  # metadata adjusting
  mt <- sapply(colData(fcs.SE), function(x) as.numeric(as.factor(x)))
  
  # prepare dr object
  if(class(dr.object) != "list"){
    diff <- setdiff(colnames(dr.object), c(colnames(colData(fcs.SE)), rownames(fcs.SE)))
    dr <- dr.object[,colnames(dr.object) %in% diff]
  }else{
    diff <- setdiff(colnames(dr.object$dr), c(colnames(colData(fcs.SE)), rownames(fcs.SE)))
    dr <- dr.object$dr[,colnames(dr.object$dr) %in% diff]
  }
  
  # combine all together
  to_export <- cbind(t(assay(fcs.SE, i = "raw")), mt, dr)
  
  ## create FCS
  flowCore::write.FCS(premessa::as_flowFrame(as.matrix(to_export), source.frame = NULL), output.name)
}                                             
