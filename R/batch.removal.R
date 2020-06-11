#' Remove batch effect
#'
#' It removes posible batch effects present in multiple FCS stored in a \code{FCS.SCE} object. It based on: \enumerate{
#'     \item Seurat's Canonical Correlation Analysis (CCA), see the \href{https://satijalab.org/seurat/Seurat_AlignmentTutorial.html}{tutorial} for more information.
#'     \item \href{https://portals.broadinstitute.org/harmony/}{Harmony}'s methodology.}
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"transformed"}.
#' @param method Methodology to perform batch correction. Possible values are "CCA" or "harmony".
#' @param batch Variable name from \code{colData(FCS.SCE)} object to remove batch effect.
#' @param new.matrix.name New normalized matrix name (it will stored within the \code{fcs.SCE} object). Default = \code{"normalized"}.
#' @param nfeatures Number of features (markers) to select as top variable features. Default = \code{10}.
#' @param dims Number of dimensions to calculated with CCA's method. Default = \code{1:10}.
#' @param threads Number of threads for multithreading. Default = \code{NULL} (i.e., all cores minus one).
#' @keywords dataset alignment normalization
#' @keywords canonical correlation
#' @keywords batch effect removal
#' @export batch.removal
#' @import future
#' @import Seurat
#' @importFrom harmony HarmonyMatrix
#' @importFrom SummarizedExperiment assay assayNames colData
#' @examples
#' \dontrun{
#'  	fcs <- batch.removal(fcs, method = "harmony", batch = "patient_id", new.matrix.name = "harmony")
#'  	fcs <- batch.removal(fcs, method = "cca", batch = "patient_id", new.matrix.name = "cca")
#' }

batch.removal <- function(fcs.SCE, assay.i = "transformed", method, batch, new.matrix.name = "normalized", nfeatures = 10, dims = 1:10, threads = NULL){
  if(tolower(method) == "cca"){
    ## setting multithreading
    library(future)
    if(is.null(threads)) threads <- availableCores()-1
    plan("multiprocess", workers = threads); options(future.globals.maxSize = 1000 * 1024^2) #https://satijalab.org/seurat/v3.0/future_vignette.html

    data <- as.Seurat(fcs.SCE, counts = assayNames(fcs.SCE)[1], data = assay.i)
    batches <- SplitObject(data, split.by = batch)
    batches <- lapply(batches, function(x) FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE))

    features <- SelectIntegrationFeatures(object.list = batches)
    batches <- lapply(batches, function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
    })

    anchors <- FindIntegrationAnchors(object.list = batches, dims = dims, reduction = "rpca")
    seurat_intg <- IntegrateData(anchorset = anchors, dims = dims)

    seurat_intg2 <- as.matrix(seurat_intg@assays$integrated@data)
    rownames(seurat_intg2) <- gsub("-", "_", rownames(seurat_intg2))
    assay(fcs.SCE, i = new.matrix.name) <- seurat_intg2[match(rownames(fcs.SCE), rownames(seurat_intg2)),match(colnames(fcs.SCE), colnames(seurat_intg2))]
    return(fcs.SCE)
  }else if(tolower(method) == "harmony"){
    norm_data <- HarmonyMatrix(data_mat = assay(fcs.SCE, assay.i), meta_data = colData(fcs.SCE),
                               vars_use = batch, do_pca = F, verbose = F)
    assay(fcs.SCE, i = new.matrix.name) <- norm_data
    return(fcs.SCE)
  }else{
    stop("Please, indicate a valid method: CCA or harmony.", call. = F)
  }
}
