#' sub.samples.idx
#'
#' This function generates an index taking (in a random manner) a number of given events by sample.
#' @param data A data.frame with the marker expression values (and metadata information binded) for each cell.
#' @param colname_samples Colname of the column with samples indentifiers.
#' @param samples_names A vector with sample identifiers.
#' @param subsampling Number of events to keep. Default = 1000
#' @param set.seed \code{\link[base:set.seed]{set.seed()}} for random picking. Default = 333
#' @param verbose Logical indicating if display process output. Default = \code{TRUE}
#' @keywords subsampling
#' @export
#' @examples
#' \dontrun{
#' sub_idx_som <- sub.samples.idx(data = mdsc_som, colname_samples = "sample_id",
#'     samples_names = md$sample_id, subsampling = 1000)
#' sel_expr <- mdsc_som[sub_idx_som,]}

sub.samples.idx <- function(data, colname_samples, samples_names, subsampling = 1000, set.seed = 333,
                            verbose = TRUE){
  set.seed(set.seed)

  sub_idx <- vector()
  for(sam in samples_names){

    if(verbose){
      print(paste0("Extracting subsampling index for: ", sam))
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
    print("Markers from panel are equal to the read FCS!")
  }else{
    print(paste0("The following marker has a different name or is not in the FCS file: ", (setdiff(markers1, markers2))))
  }
}


#' lims.FCS
#'
#' This function estimates the minimum and maximum values for an object with multiple FCS files. Usefull for plotting limits.
#' @param fcs A data.frame with marker expression values for each cell.
#' @export
#' @examples
#' \dontrun{lims.FCS(fcs_raw)}

lims.FCS <- function(fcs){
  mn <- c()
  mx <- c()

  for(i in names(fcs@frames)){
    aux <- as.data.frame(summary(flowCore::exprs(fcs@frames[[i]])))
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
#' @examples
#' \dontrun{
#' markers.names(fcs_raw)
#' fcs_raw <- markers.names(fcs_raw, new_names = c("FSC_A", "FSC_H", "SSC_A", "SSC_H",
#'     "CD62L", "CXCR3", "CD8", "CCR4", "CCR6", "CD4", "CD27", "CD45RA"))}


markers.names <- function(fcs, new_names = NULL){
  if(is.null(new_names)){
    gridExtra::grid.arrange(gridExtra::tableGrob(data.frame(fcs[[1]]@parameters@data$name, fcs[[1]]@parameters@data$desc, "not defined yet!"),
                           cols = c("Florophore/Channel", "Marker (default name)", "New marker name")))
  }else{
    gridExtra::grid.arrange(gridExtra::tableGrob(data.frame(fcs[[1]]@parameters@data$name, fcs[[1]]@parameters@data$desc, new_names),
                           cols = c("Florophore/Channel", "Marker (default name)", "New marker name")))
    colnames(fcs) <- new_names
    # suppressWarnings(assign(deparse(substitute(fcs)), aux_fcs))
  }
  return(fcs)
}
