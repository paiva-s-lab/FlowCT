#' Show and change markers names
#'
#' It shows marker names of a \code{fcs.SCE} object and renames them accordin a new vector provided by the user.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param new.names Vector with new channel/marker names (it must has the same length that \code{fcs.SCE}'s markers). Default = \code{NULL} (i.e., markers will not be renamed, only displayed).
#' @keywords marker renaming names
#' @importFrom knitr kable
#' @export marker.names
#' @examples
#' \dontrun{
#' marker.names(fcs)
#' new_names <- c("FSC_A", "FSC_H", "SSC_A", "SSC_H", "CD62L", 
#'      "CXCR3", "CD8", "CCR4", "CCR6", "CD4", "CD45", "CD27")
#' fcs <- marker.names(fcs, new.names = new_names)
#' }

marker.names <- function (fcs.SCE, new.names = NULL){
  raw_names <- rownames(fcs.SCE)
  if(is.null(new.names)){
    if(is.null(fcs.SCE@metadata$raw_channel_names)){
      df <- data.frame(raw_names, "not defined yet!")
      colnames(df) <- c("raw name", "new name")
      asciify(df)
    }else{
      df <- data.frame(fcs.SCE@metadata$raw_channel_names, raw_names)
      colnames(df) <- c("raw name", "new name")
      asciify(df)
    }
  }else{
    ## add new markers names
    if(length(new.names) != length(raw_names)) stop("New names must to have the same length to the original ones!")
    df <- data.frame(raw_names, new.names)
    colnames(df) <- c("raw name", "new name")
    asciify(df)
    
    rownames(fcs.SCE) <- new.names
    fcs.SCE@metadata$raw_channel_names <- raw_names
    return(fcs.SCE)
  }
}
