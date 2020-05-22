#' sub.samples
#'
#' It generates reduce randomly the number of events of a fcs.SCE object (it computes this reduction for each FCS file separatelly inside this object). It can generate a new reduced \code{fcs.SCE} object or a simple index position for removing events.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{fcs.SCE()}}.
#' @param subsampling Number of events to keep. Default = 1000.
#' @param index Logical indicating if returning an fcs.SCE object or a index vector. Default = \code{FALSE}.
#' @keywords subsampling
#' @export
#' @importFrom progress progress_bar
#' @examples
#' \dontrun{
#' # option 1
#' fcs_red <- sub.samples(data = fcs, subsampling = 1000)
#'
#' # option 2
#' keep_red <- sub.samples(data = fcs, subsampling = 1000, index = T)
#' fcs_red <- fcs[,keep_red]
#' }

sub.samples <- function (fcs.SCE, subsampling = 1000, index = F){
  set.seed(333)
  pb <- progress_bar$new(total = length(unique(fcs.SCE$filename)), format = "Random subsampling [:bar]")
  
  sub_idx <- vector()
  for(i in unique(fcs.SCE$filename)) {
    pb$tick()
    
    aux <- colnames(fcs.SCE[,fcs.SCE$filename == i])
    if(length(aux) < subsampling){
      cat("Attention!: filename", i, "has a lower number of events, reduction won't be computed.\n")
      sub_idx <- append(sub_idx, aux)
    }else{
      sub_idx <- append(sub_idx, aux[sample(length(aux), subsampling)])    
    }
    
    Sys.sleep(1/10)
  }
  if(index) return(sub_idx) else return(fcs.SCE[,sub_idx])
}
