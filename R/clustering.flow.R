#' Clustering FCS data
#'
#' This function offers multiple automatic clustering methods for a \code{FCS.SCE} object. Those methods available are:\enumerate{
#'           \item SOM: Self-Organizing Map (SOM) based on \href{https://github.com/SofieVG/FlowSOM}{\code{FlowSOM} package}.
#'           \item Phenograph: An implementation of PhenoGraph algorithm (\href{https://github.com/JinmiaoChenLab/Rphenograph}{here} for more information).
#'           \item Seurat: Clustering method based on \href{https://satijalab.org/seurat/}{Seurat}'s single-cell analysis.
#'           \item PARC: R implementation for Phenotyping by Accelerated Refined Community-partitioning (PARC) method, see \href{https://github.com/ShobiStassen/PARC}{Python module}.}
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate SOM clustering. Default = \code{"normalized"}.
#' @param method What method should be used for clustering purposes. Available ones are "SOM", "Phenograph", "Seurat" and "PARC".
#' @param scale Should be data be scale before SOM clustering? (only available for this method). Default = \code{"FALSE"}.
#' @param markers.to.use Markers used for considering within the clustering calculation. Default = \code{"all"}.
#' @param num.k Number of clusters to calculate for methods "SOM" and "Phenograph".
#' @param seurat.res Seurat's resolution to calculate clustering (it indicates the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters). Default = \code{0.4}.
#' @param seurat.dims Number of dimensions to calculated with Seurat's method. Default = \code{1:10}.
#' @keywords SOM Seurat PARC Phenograph unsupervised clustering
#' @export clustering.flow
#' @import FlowSOM
#' @import Seurat
#' @import reticulate
#' @importFrom Rphenograph Rphenograph
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @importFrom SummarizedExperiment colData assay
#' @examples
#' \dontrun{
#' fcs <- clustering.flow(fcs, method = "som", num.k = 20)
#' fcs <- clustering.flow(fcs, method = "phenograph", num.k = 40)
#' fcs <- clustering.flow(fcs, method = "Seurat", seurat.res = 0.5)
#' fcs <- clustering.flow(fcs, method = "parc", assay.i = "transformed")
#' }


clustering.flow <- function(fcs.SCE, assay.i = "normalized", method, scale = FALSE, markers.to.use = "all", num.k, seurat.res = 0.4, seurat.dims = 1:10){
  set.seed(333)

  if(tolower(method) == "som"){
    data <- as.flowSet.SE(fcs.SCE, assay.i)
    if(length(markers.to.use) == 1 && markers.to.use == "all") markers.to.use <- rownames(data)

    fsom <- suppressMessages(ReadInput(data, transform = FALSE, scale = scale)) #read data
    fsom <- suppressMessages(BuildSOM(fsom, colsToUse = markers.to.use)) #build SOM

    ## metaclustering
    mc <- suppressMessages(ConsensusClusterPlus(t(fsom$map$codes), maxK = num.k, reps = 100,
                                                pItem = 0.9, pFeature = 1, title = "consensus_plots", plot = "pdf",
                                                clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                                                distance = "euclidean", seed = 333, verbose = F))
    unlink("consensus_plots", recursive = TRUE)

    code_clustering1 <- mc[[num.k]]$consensusClass %>% as.factor() #get cluster ids for each cell
    cell_clustering1 <- code_clustering1[fsom$map$mapping[,1]]

    colData(fcs.SCE)[,paste0("SOM.k", num.k)] <- cell_clustering1
    return(fcs.SCE)
  }else if(tolower(method) == "seurat"){
    data <- as.Seurat(fcs.SCE, counts = assayNames(fcs.SCE)[1], data = assay.i)
    if(length(markers.to.use) == 1 && markers.to.use == "all") markers.to.use <- rownames(data)

    s <- ScaleData(data, verbose = FALSE)
    s <- RunPCA(s, features = markers.to.use, verbose = FALSE)
    s <- FindNeighbors(s, dims = seurat.dims)
    s <- FindClusters(s, resolution = seurat.res)

    colData(fcs.SCE)[,paste0("seurat.r", seurat.res)] <- as.factor(as.numeric(as.character(s$seurat_clusters))+1)
    return(fcs.SCE)
  }else if(tolower(method) == "phenograph"){
    data <- as.matrix(t(assay(fcs.SCE, assay.i)))
    if(length(markers.to.use) == 1 && markers.to.use == "all") markers.to.use <- colnames(data)

    Rphenograph_out <- Rphenograph(data[,markers.to.use], k = num.k)
    colData(fcs.SCE)[,paste0("phenog.k", num.k)] <- factor(Rphenograph_out[[2]]$membership)
    return(fcs.SCE)
  }else if(tolower(method) == "parc"){
    p <- import("parc")
    pres <- p$PARC(t(assay(fcs.SCE, assay.i)))
    pres$run_PARC()
    fcs.SCE$PARC <- as.factor(as.numeric(as.character(unlist(pres$labels)))+1)
    return(fcs.SCE)
  }else{
    stop("Please, indicate a valid method: SOM, Seurat, Phenograph or PARC", call. = F)
  }
}
