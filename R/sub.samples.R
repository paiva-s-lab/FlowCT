#' Reduce the nummber of events of a \code{fcs.SCE} object
#'
#' It generates reduce randomly the number of events of a fcs.SCE object (it computes this reduction for each FCS file separatelly inside this object). It can generate a new reduced \code{fcs.SCE} object or a simple index position for removing events.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param subsampling Number of events to keep (in each FCS file contained within the \code{fcs.SCE} object). If indicated value is between 0 and 1, it will considered as a percentage to keep instead an absolute number. Default = 1000.
#' @param index Logical indicating if returning an fcs.SCE object or a index vector. Default = \code{FALSE}.
#' @param unbalanced If indicated, subsampling will be applied in a "unbalanced" way, i.e., only in the specified element in a given column. Its format is a vector containing the column name to reduce and the element within that column to be subsampled. It is very IMPORTANT keep in mind this subsampling is only for visualization purposes, not for downstream statistical analysis. Default = \code{NULL}.
#' @keywords subsampling
#' @export
#' @importFrom progress progress_bar
#' @examples
#' \dontrun{
#' # option 1, simple subsampling
#' fcs_red <- sub.samples(data = fcs, subsampling = 1000)
#'
#' # option 2, return the subsampling index
#' keep_red <- sub.samples(data = fcs, subsampling = 1000, index = T)
#' fcs_red <- fcs[,keep_red]
#'
#' # option 3, unbalanced subsampling (in percentaje way)
#' fcs_red <- sub.samples(data = fcs, subsampling = 0.3, unbalanced = c("condition", "not_progressed")
#' }

sub.samples <- function (fcs.SCE, subsampling = 1000, index = F, unbalanced = NULL){
  set.seed(333)
  pb <- progress_bar$new(total = length(unique(fcs.SCE$filename)), format = "Random subsampling [:bar]")
  sub_idx <- vector()

  if(!is.null(unbalanced)){
    aux0 <- fcs.SCE[,fcs.SCE[[unbalanced[1]]] == unbalanced[2]]
    aux_add <- fcs.SCE[,fcs.SCE[[unbalanced[1]]] != unbalanced[2]]

    for(i in unique(aux0$filename)){
      pb$tick()
      
      aux <- colnames(aux0[,aux0$filename == i])
      if(is.null(aux)) {stop("fcs.SCE object is not correctly generated, missing single-cell indentifiers...")}
      if(between(subsampling, 0, 1)) subsampling0 <- length(aux)*subsampling
      
      if(length(aux) < subsampling0){
        cat("Attention!: file", i, "has a lower number of events, reduction won't be computed.\n")
        sub_idx <- append(sub_idx, aux)
      }else{
        sub_idx <- append(sub_idx, aux[sample(length(aux), subsampling0)])    
      }
        Sys.sleep(1/10)
    }
    cat("\n")

    if(index) return(sub_idx) else return(cbind(aux0[,sub_idx], aux_add))

  }else{
    for(i in unique(fcs.SCE$filename)) {
      pb$tick()
      
      aux <- colnames(fcs.SCE[,fcs.SCE$filename == i])
      if(is.null(aux)) {stop("fcs.SCE object is not correctly generated, missing single-cell indentifiers...")}
      if(between(subsampling, 0, 1)) subsampling0 <- length(aux)*subsampling
      
      if(length(aux) < subsampling0){
        cat("Attention!: file", i, "has a lower number of events, reduction won't be computed.\n")
        sub_idx <- append(sub_idx, aux)
      }else{
        sub_idx <- append(sub_idx, aux[sample(length(aux), subsampling0)])    
      }
      Sys.sleep(1/1000)
    }
   cat("\n")

   if(index) return(sub_idx) else return(fcs.SCE[,sub_idx]) 
  }
}
