#' sub.samples
#'
#' It generates reduce randomly the number of events of a FCS.SE object (it computes this reduction for each FCS file separatelly inside this object). It can generate a new reduced FCS.SE object or a simple index position for removing events.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param subsampling Number of events to keep. Default = 1000.
#' @param index Logical indicating if returning an FCS.SE object or a index vector. Default = \code{FALSE}.
#' @keywords subsampling
#' @export
#' @examples
#' \dontrun{
#' # option 1
#' fcs_se_red <- sub.samples(data = fcs_se, subsampling = 1000)
#'
#' # option 2
#' keep_red <- sub.samples(data = fcs_se, subsampling = 1000, index = T)
#' fcs_se_red <- fcs_se[,keep_red]
#' }

sub.samples <- function (fcs.SE, subsampling = 1000, index = F){
  set.seed(333)
  pb <- progress::progress_bar$new(total = length(unique(fcs.SE$filename)), format = "Random subsampling [:bar]")
  
  sub_idx <- vector()
  for(i in unique(fcs.SE$filename)) {
    pb$tick()
    
    aux <- colnames(fcs.SE[,fcs.SE$filename == i])
    if(length(aux) < subsampling){
      cat("Attention!: filename", i, "has a lower number of events, reduction won't be computed.\n")
      sub_idx <- append(sub_idx, aux)
    }else{
      sub_idx <- append(sub_idx, aux[sample(length(aux), subsampling)])    
    }
    
    Sys.sleep(1/10)
  }
  if(index) return(sub_idx) else return(fcs.SE[,sub_idx])
}
