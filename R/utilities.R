#' sub.samples.idx
#'
#' This function generates an index taking (in a random manner) a number of given events by sample. Later, this index can be used to subsampling an object (ie. a data.frame)
#' @param data A data.frame with the marker expression values (and metadata information merged) for each cell.
#' @param colname_samples Name of the column with sample identifiers to generate the index.
#' @param subsampling Number of events to keep. Default = 1000
#' @param set.seed \code{\link[base:set.seed]{set.seed()}} for random picking. Default = 333
#' @param verbose Logical indicating if display process output. Default = \code{TRUE}
#' @keywords subsampling
#' @export
#' @examples
#' \dontrun{
#' sub_idx_som <- sub.samples.idx(data = mdsc_som, colname_samples = "sample_id",
#'     subsampling = 1000)
#' sel_expr <- mdsc_som[sub_idx_som,]}

sub.samples.idx <- function(data, colname_samples, subsampling = 1000, set.seed = 333,
                            verbose = TRUE){
  set.seed(set.seed)
  
  sub_idx <- vector()
  for(sam in unique(colname_samples)){
    
    if(verbose){
      cat("Extracting subsampling index for: ", sam, "\n", sep = "")
    }
    
    aux <- rownames(data[data[,colname_samples] == sam,])
    sub_idx <- append(sub_idx, aux[sample(length(aux), subsampling)])
  }
  
  sub_idx <- as.integer(sub_idx)
  return(sub_idx)
}


#' markers.equal
#'
#' This function compares if two strings are equal (in this case, for marker comparison).
#' @param markers1 First string with the marker names.
#' @param markers2 Second string with the marker names.
#' @export
#' @examples
#' \dontrun{markers.equal(surface_markers, fcs_raw@colnames)}

markers.equal<- function(markers1, markers2){
  if(all(markers1 %in% markers2) == TRUE){
    cat("Markers from panel are equal to the read FCS!\n")
  }else{
    cat("The following marker has a different name or is not in the FCS file: ", (setdiff(markers1, markers2)), "\n", sep = "")
  }
}



#' lims.FCS
#'
#' This function estimates the minimum and maximum values for an object with multiple FCS files. Useful for plotting limits.
#' @param fcs A \code{flowSet} object read with \code{\link[flowCore:read.flowSet]{flowCore::read.flowSet()}}.
#' @export
#' @importFrom flowCore exprs
#' @examples
#' \dontrun{lims.FCS(fcs_raw)}

lims.FCS <- function(fcs){
  mn <- mx <- c()
  
  for(i in names(fcs@frames)){
    aux <- as.data.frame(summary(exprs(fcs@frames[[i]])))
    aux2 <- data.frame(marker = aux$Var2,
                       statistic = sapply(aux$Freq, function(x) strsplit(as.character(x), ":")[[1]][1]),
                       value = as.numeric(sapply(aux$Freq, function(x) strsplit(as.character(x), ":")[[1]][2])))
    aux_mx <- max(as.numeric(aux2$value[grep("Max", aux2$statistic)]))
    aux_mn <- min(as.numeric(aux2$value[grep("Min", aux2$statistic)]))
  }
  mn <- append(mn, aux_mn)
  mx <- append(mx, aux_mx)
  
  return(c(mn, mx))
}


#' markers.names
#'
#' This function shows default names for chanels in FCS files and offers the option to change them manually.
#' @param fcs A \code{flowSet} object read with \code{\link[flowCore:read.flowSet]{flowCore::read.flowSet()}}.
#' @param new_names A string with new markers names to rename the \code{flowSet} object. Default = \code{NULL}
#' @export
#' @import gridExtra
#' @importFrom ggplot2 ggsave
#' @importFrom grid grid.draw
#' @examples
#' \dontrun{
#' markers.names(fcs_raw) #it only shows marker names...
#' fcs_raw <- markers.names(fcs_raw, new_names = c("FSC_A", "FSC_H", "SSC_A", "SSC_H",
#'     "CD62L", "CXCR3", "CD8", "CCR4", "CCR6", "CD4", "CD27", "CD45RA"))} #it shows and changes marker names according supplied ones

markers.names <- function (fcs, new_names = NULL){
  # require(gridExtra)
  lapply(c("grid", "gridExtra", "ggplot2"), require, character.only = TRUE)

  if(is.null(new_names)){
    if(all(fcs[[1]]@parameters@data$name == fcs@colnames)){
      print(knitr::kable(data.frame(fcs[[1]]@parameters@data$name, 
                                                  fcs[[1]]@parameters@data$desc, "not defined yet!"), 
                                       col.names = c("Florophore/Channel", "Marker (default name)", 
                                                "New marker name")))
    }else{
      print(knitr::kable(data.frame(fcs[[1]]@parameters@data$name, 
                                  fcs[[1]]@parameters@data$desc, fcs@colnames), 
                                       col.names = c("Florophore/Channel", "Marker (default name)", "New marker name")))
    }

  ## add new markers names
  }else{
    print(knitr::kable(data.frame(fcs[[1]]@parameters@data$name, 
                                                fcs[[1]]@parameters@data$desc, new_names), 
                                     col.names = c("Florophore/Channel", "Marker (default name)", "New marker name")))

    fcs@colnames <- new_names
    return(fcs)
  }
}
