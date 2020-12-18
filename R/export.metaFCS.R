#' Create FCS from a \code{fcs.SCE} object
#'
#' It creates a FCS file containing all analysis incorpored to colData(fcs.SCE) as well as the dimensional reduction coordinates calculated.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param dr.object Object created with \code{\link[FlowCT:dim.reduction]{FlowCT::dim.reduction()}} function or a table combining DR and experimental metadata information.
#' @param output.name Name for generated FCS file. If \code{NULL} (default), the project name specified in the \code{fcs.SCE} building will be used.
#' @param output.folder Name of the folder within store the new generated FCS files. Default = \code{getwd()} (i.e., current directory).
#' @param separate.fcs Logical indicating whether FCS must be exporting separatelly or not. Default = \code{FALSE}.
#' @param output.suffix Suffix to be added to the new generated FCS files. Default = \code{"meta"}.
#' @keywords FCS generation
#' @keywords final FCS
#' @keywords exporting
#' @export
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay
#' @importFrom premessa as_flowFrame
#' @importFrom flowCore write.FCS
#' @examples
#' \dontrun{
#' # option 1: save entire analysis as a single FCS file
#' export.metaFCS(fcs.SCE = fcs, output.name = "final_analysis")
#'
#' # option 2: save separatelly each file contained within the fcs.SCE object
#' export.metaFCS(fcs.SCE = fcs, separate.fcs = T, output.suffix = "meta.final")
#' }

export.metaFCS <- function(fcs.SCE, assay.i = "raw", output.name = NULL, output.folder = getwd(), separate.fcs = FALSE, output.suffix = "meta"){
  if(!separate.fcs){
    # metadata adjusting
    mt <- sapply(colData(fcs.SCE), function(x) as.numeric(as.factor(x)))

    # prepare dr object
    if(length(fcs.SCE@int_colData@listData$reducedDims) > 0){
      drs <- as.data.frame(fcs.SCE@int_colData@listData$reducedDims@listData)
    }else{
      drs <- NULL
    }

    # combine all together
    to_export <- cbind(t(assay(fcs.SCE, i = assay.i)), mt, drs)

    ## create FCS
    if(is.null(output.name)) output.name <- fcs.SCE@metadata$project_name
    extension <- tail(strsplit(fcs.SCE@metadata$input_fcs[1], "\\.")[[1]], 1)

    write.FCS(as_flowFrame(as.matrix(to_export)), paste0(output.folder, "/", output.name, ".", extension))
  }else{
    for(i in unique(fcs.SCE$filename)){
      aux_fcs.SCE <- fcs.SCE[,fcs.SCE$filename == i]
      mt <- sapply(colData(aux_fcs.SCE), function(x) as.numeric(as.factor(x)))

      if(length(aux_fcs.SCE@int_colData@listData$reducedDims) > 0){
        drs <- as.data.frame(aux_fcs.SCE@int_colData@listData$reducedDims@listData)
      }else{
        drs <- NULL
      }

      to_export <- cbind(t(assay(aux_fcs.SCE, i = assay.i)), mt, drs)
      extension <- tail(strsplit(basename(i), "\\.")[[1]], 1)
      filename <- gsub(paste0(".", extension), "", basename(i))

      cat("Saving new: ", filename, ".meta.", extension, "\n", sep = "")
      write.FCS(as_flowFrame(as.matrix(to_export)), paste0(output.folder, "/", filename, ".", output.suffix, ".", extension))
    }
  }
}
