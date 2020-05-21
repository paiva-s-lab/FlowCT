#' fcs.SCE
#'
#' It reads and creates a fcs.SCE object (based on the structure of a (\ref{https://www.bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html}{\code{SummarizedExperiment}}) from FCS files in a specific folder or indicated in a vector.
#' @param filelist A vector with full path of FCS files to be read, commonly generated through \code{\link[base:list.files]{base::list.files()}}. If \code{NULL}, this file list will be generated as indicated below.
#' @param directory If \code{filelist = NULL}, those files stored in this location will be read. Default = \code{getwd()} (current directory).
#' @param pattern Pattern for reading files within \code{directory}. Default = \code{"fcs"}.
#' @param events Numeric vector indicating how many events to read in each FCS file (same behaviour as \code{\link[flowCore:read.FCS]{flowCore::read.FCS()}}), those are the number of events for new and unified generated FCS files. Default = \code{"all"}.
#' @param dataset Specify data segment to be read in the FCS file (same behaviour as \code{\link[flowCore:read.FCS]{flowCore::read.FCS()}}). Default = \code{1}.
#' @param num.threads Number of threads for reading the FCS files. For Windows users, you should to install \code{parallelsugar} from (\href{https:https://github.com/nathanvan/parallelsugar}{Nathanvan's GitHub}). Defult = \code{NULL}, i.e., not parallelization.
#' @param metadata Metadata information for each FCS files included in the analysis. It will be used to assign to each event an specific metadata and it will be added to \code{fcs.SCE@metadata} as \code{reduced_metadata}.
#' @param transformation Because FCS data is normally exported in a logarithmic format, it is necessary its transformation for downstream analysis. Default (and the only one) = \code{arcsinh}.
#' @param transf.cofactor Cofactor numeric value for that \code{arcsinh} transformation. Values can ranged from 15 to 50 (for mass cytometry) and from 150 to 10,000 (for flow cytometry). Default = \code{500}.
#' @param project.name String naming the \code{fcs.SCE} generated. Default = \code{"noname"}.
#' @keywords fcs.SCE object
#' @keywords read FCS
#' @keywords single-cell metadata
#' @return The fcs.SCE object is comprised with multiple (and downstream customizable) elements \enumerate{
#'           \item \code{colData(fcs.SCE)}: A \code{data.frame} with metadata information for each event in the read FCS files.
#'           \item \code{assayNames(fcs.SCE)}: Two expression matrices with markers as rows and events (cells) as columns. The first matrix (\code{raw}) contains the raw expression for each marker, and the second one (\code{transformed}), the \code{arcsinh}-transformed expression for subsequent analysis.
#'           \item \code{metadata(fcs.SCE)}: Additional information about the experiment and following downstream steps: those read FCS files (\code{input_files}), \code{reduced_metada} with metadata for each FCS file...
#' @export
#' @importFrom flowCore fsApply exprs
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @examples
#' \dontrun{
#' head(metadata_user)
#'  #       filename sample_id condition patient_id
#'  # 1 BM_017294.fcs         1        BM     017294
#'  # 2 BM_017564.fcs         2        BM     017564
#'  # 3 BM_017612.fcs         3        BM     017612
#'  # 4 PB_017294.fcs         4        PB     017294
#'  # 5 PB_017564.fcs         5        PB     017564
#'  # 6 PB_017612.fcs         6        PB     017612
#' fcs <- fcs.SCE(directory = "../data/", which.lines = 1000, metadata = metadata_user, transf.cofactor = 500)
#' }

fcs.SCE <- function(filelist = NULL, directory = getwd(), pattern = ".fcs$", events = "all", dataset = 1, num.threads = NULL, metadata, transformation = "arcsinh", transf.cofactor = 500, project.name = "noname"){
  # read FCS files
  fcs <- suppressMessages(read.FCSset(filelist, directory, pattern, events = events, dataset, num.threads))
  raw_data <- fsApply(fcs, exprs)
  colnames(raw_data) <- paste0(fcs[[1]]@parameters@data$name, ":", fcs[[1]]@parameters@data$desc)
  
  # build sc metada
  col_data <- suppressMessages(scMetadata.fromFCS(fcs, metadata, add.exprs = F))
  
  # transformation
  if(!is.null(transformation)){
    if(transformation != "arcsinh"){
      cat("Sorry, the only transformation method available (yet) is 'arcsinh'")
    }
    transf_data <- asinh(raw_data/transf.cofactor)
    colnames(transf_data) <- paste0(fcs[[1]]@parameters@data$name, ":", fcs[[1]]@parameters@data$desc)      
  }
  
  # build SingleCellExperiment object
  fcs.SCE <- SingleCellExperiment(assays = list(raw = t(raw_data), transformed = t(transf_data)), colData = col_data, 
                          # additional info
                          metadata = list(project_name = project.name, 
                                          input_fcs = as.vector(unique(fcs@phenoData@data$name)), 
                                          reduced_metadata = metadata))
 
  return(fcs.SCE)  
}
