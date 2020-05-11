# 'marker.names
#'
#' It shows marker names of a \code{FCS.SE} object and renames them accordin a new vector provided by the user.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param new.names Vector with new channel/marker names (it must has the same length that \code{FCS.SE}'s markers). Default = \code{NULL} (i.e., markers will not be renamed, only displayec).
#' @keywords median values
#' @keywords MFI
#' @keywords median fluorescence intensity
#' @export median.values
#' @examples
#' \dontrun{
#' marker.names(fcs_se)
#' new_names <- c("FSC_A", "FSC_H", "SSC_A", "SSC_H", "CD62L", "CXCR3", "CD8", "CCR4", "CCR6", "CD4", "CD45", "CD27")
#' fcs_se <- marker.names(fcs_se, new.names = new_names)
#' }

marker.names <- function (fcs.SE, new.names = NULL){
  raw_names <- rownames(fcs.SE)
  if(is.null(new.names)){
    if(is.null(fcs.SE@metadata$raw_channel_names)){
      print(knitr::kable(data.frame(raw_names, "not defined yet!"), col.names = c("raw name", "new name")))
    }else{
      print(knitr::kable(data.frame(fcs.SE@metadata$raw_channel_names, raw_names), 
                         col.names = c("raw name", "new name")))
    }
  }else{
    ## add new markers names
    if(length(new.names) != length(raw_names)) stop("New names must to have the same length to the original ones!")
    print(knitr::kable(data.frame(raw_names, new.names), col.names = c("raw name", "new name")))
    
    rownames(fcs.SE) <- new.names
    fcs.SE@metadata$raw_channel_names <- raw_names
    return(fcs.SE)
  }
}
