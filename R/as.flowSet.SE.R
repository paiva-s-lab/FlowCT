#' Transform a \code{fcs.SCE} object into a \code{flowset} object
#'
#' It tranforms a \code{fcs.SCE} object into a \href{https://www.rdocumentation.org/packages/flowCore/versions/1.38.2/topics/flowSet-class}{\code{flowset}} object.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which obtain a \code{flowset} object.
#' @keywords fcs.SCE to flowset
#' @keywords flowset generation
#' @export
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay
#' @importFrom premessa as_flowFrame
#' @importFrom methods as

as.flowSet.SE <- function(fcs.SCE, assay.i){
  df_flowset <- lapply(unique(fcs.SCE$filename), function(x) as_flowFrame.me(t(assay(fcs.SCE[,fcs.SCE$filename == x], i = assay.i))))
  names(df_flowset) <- unique(fcs.SCE$filename)
  return(as(df_flowset, "flowSet"))
}

# modified from https://github.com/ParkerICI/premessa/blob/master/R/fcs_io.R
as_flowFrame.me <- function(exprs.m) {
    flow.frame <- flowCore::flowFrame(exprs.m)
	params <- flowCore::parameters(flow.frame)
    pdata <- flowCore::pData(params)

    for (i in 1:ncol(flow.frame)) {
        s <- paste("$P",i,"S",sep="")
        n <- paste("$P",i,"N",sep="")
        r <- paste("$P",i,"R",sep="")
        b <- paste("$P",i,"B",sep="")
        e <-  paste("$P",i,"E",sep="")

        keyval <- list()

        keyval[[s]] <- colnames(exprs.m)
        keyval[[n]] <- colnames(exprs.m)[i]
        keyval[[r]] <- ceiling(max(exprs.m[,i], na.rm = TRUE))

        keyval[[b]] <- 32
        keyval[[e]] <- "0,0"
        flowCore::keyword(flow.frame) <- keyval

        pdata[i,"minRange"] <- min(exprs.m[,i], na.rm = TRUE)
        pdata[i,"maxRange"] <- max(exprs.m[,i], na.rm = TRUE)
    }

    flowCore::pData(params) <- pdata
    flowCore::parameters(flow.frame) <- params

    return(flow.frame)
}
