#' Remove batch effect
#'
#' It removes posible batch effects present in multiple FCS stored in a \code{FCS.SCE} object. It based on: \enumerate{
#'     \item Seurat's Canonical Correlation Analysis (CCA) or Reciprocal PCA (RPCA), see the \href{https://satijalab.org/seurat/Seurat_AlignmentTutorial.html}{tutorial} for more information (another \href{https://satijalab.org/seurat/v3.0/integration.html}{one}).
#'     \item \href{https://portals.broadinstitute.org/harmony/}{Harmony}'s methodology.}
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"transformed"}.
#' @param method Methodology to perform batch correction. Possible values are "Seurat" or "harmony".
#' @param batch Variable name from \code{colData(FCS.SCE)} object to remove batch effect.
#' @param new.matrix.name New normalized matrix name (it will stored within the \code{fcs.SCE} object). Default = \code{"normalized"}.
#' @param harmony.params List with those options from original \code{\link[harmony:HarmonyMatrix]{harmony::HarmonyMatrix()}} function that user wants to customize. Default, harmony's default ones.
#' @param seurat.params List with those options from original \code{\link[Seurat:FindIntegrationAnchors]{FindIntegrationAnchors()}} function that user wants to customize. Default, Seurat's default ones except \code{dims = 1:10} and \code{reduction = "rpca"} (the other option is "cca", but slower for typical flow cytometry dataset sizes).
#' @param threads Number of threads for multithreading. Default = \code{NULL} (i.e., all available cores minus one).
#' @keywords dataset alignment normalization
#' @keywords canonical correlation
#' @keywords batch effect removal
#' @export batch.removal
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay
#' @examples
#' \dontrun{
#'  fcs <- batch.removal(fcs, method = "harmony", batch = "patient_id", 
#'    new.matrix.name = "harmony")
#'  fcs <- batch.removal(fcs, method = "harmony", batch = "patient_id", 
#'    new.matrix.name = "harmony",
#'    harmony.params = list(plot_convergence = T)) #display optimal number of iterations
#'  fcs <- batch.removal(fcs, method = "seurat", batch = "patient_id", 
#'    seurat.params = list(reduction = "cca"))
#'  fcs <- batch.removal(fcs, method = "seurat", batch = "patient_id",
#'  	              seurat.params = list(reduction = "cca", dims = 1:50))
#' }

batch.removal <- function(fcs.SCE, assay.i = "transformed", method, batch, new.matrix.name = "normalized", harmony.params = NULL, seurat.params = NULL, threads = NULL){
  if(tolower(method) == "seurat"){
    if (!requireNamespace(c("future", Seurat), quietly = TRUE)) stop("Packages \"future\" and \"Seurat\" needed for this function to work. Please install them.", call. = FALSE)
    # require(future)
    # require(Seurat)

    ## setting multithreading
    if(is.null(threads)) threads <- availableCores()-1
    plan("multiprocess", workers = threads); options(future.globals.maxSize = 1000 * 1024^2) #https://satijalab.org/seurat/v3.0/future_vignette.html

    ## prepare internal options
    seurat.defaults <- list(dims = 1:10, reduction = "rpca", k.anchor = 5, k.filter = 200,
                            k.score = 20, max.features = 200, nn.method = "rann", eps = 0)
    if(!is.null(seurat.params)){
      diffs <- setdiff(names(seurat.defaults), names(seurat.params))
      seurat.params2 <- c(seurat.params, seurat.defaults[diffs])
    }else{
      seurat.params2 <- seurat.defaults
    }

    ## processing...
    data <- as.Seurat(fcs.SCE, counts = assayNames(fcs.SCE)[1], data = assay.i)
    batches <- SplitObject(data, split.by = batch)
    batches <- lapply(batches, function(x) FindVariableFeatures(x, selection.method = "vst",
                                                                nfeatures = length(fcs), verbose = FALSE))

    features <- SelectIntegrationFeatures(object.list = batches)
    batches <- lapply(batches, function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE, npcs = length(seurat.defaults$dims))
    })


    anchors <- do.call(FindIntegrationAnchors, c(list(object.list = batches), seurat.params2))
    seurat_intg <- IntegrateData(anchorset = anchors, dims = seurat.defaults$dims)

    seurat_intg2 <- as.matrix(seurat_intg@assays$integrated@data)
    rownames(seurat_intg2) <- gsub("-", "_", rownames(seurat_intg2))
    assay(fcs.SCE, i = new.matrix.name) <- seurat_intg2[match(rownames(fcs.SCE), rownames(seurat_intg2)),
                                                        match(colnames(fcs.SCE), colnames(seurat_intg2))]
    return(fcs.SCE)
  }else if(tolower(method) == "harmony"){
    if (!requireNamespace("harmony", quietly = TRUE)) stop("Package \"harmony\" (https://github.com/immunogenomics/harmony) needed for this function to work. Please install it.", call. = FALSE)

    ## prepare internal options
    if(!is.null(harmony.params)){
      harmony.defaults <- formals(HarmonyMatrix)[-(1:4)]
      diffs <- setdiff(names(harmony.defaults), names(harmony.params))
      harmony.params2 <- c(harmony.params, harmony.defaults[diffs])
    }else{
      harmony.params2 <- harmony.defaults
    }

    norm_data <- do.call(HarmonyMatrix, c(list(data_mat = assay(fcs.SCE, assay.i), meta_data = colData(fcs.SCE),
                                               vars_use = batch, do_pca = F), harmony.params2))
    assay(fcs.SCE, i = new.matrix.name) <- norm_data
    return(fcs.SCE)
  }else{
    stop("Please, indicate a valid method: Seurat or harmony.", call. = F)
  }
}
