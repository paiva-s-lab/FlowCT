#' read.FCSset
#'
#' It reads multiple FCS files (in a computationaly-distributed way) contained in a specific folder or given through a vector.
#' @param filelist A vector with full path of FCS files to be read, commonly generated through \code{\link[base]{list.files()}}. If \code{NULL}, this file list will be generated as indicated below.
#' @param directory If \code{filelist = NULL}, those files stored in this location will be read. Default = \code{getwd()} (current directory).
#' @param pattern Pattern for reading files within \code{directory}. Default = \code{"fcs"}.
#' @param events Numeric vector indicating how many events to read in each FCS file (same behaviour as \code{\link[flowCore]{read.FCS()}}). Default = \code{NULL}, i.e., all events will be read.
#' @param dataset Specify data segment to be read in the FCS file (same behaviour as \code{\link[flowCore]{read.FCS()}}). Default = \code{1}.
#' @param num.threads Number of threads for reading the FCS files. For Windows users, you should to install \code{parallelsugar} package from (\href{https:https://github.com/nathanvan/parallelsugar}{Nathanvan's GitHub}). Defult = \code{NULL} (i.e., not parallelization).
#' @keywords FCS reading
#' @keywords FCS parallel
#' @export
#' @importFrom ncdfFlow as.flowSet read.ncdfFlowSet
#' @examples
#' \dontrun{
#' # option 1: trough a vector with filenames (full path)
#' filelist <- list.files(pattern = "fcs", path = "../data/", full.names = T)
#' fcs1 <- read.FCSset(filelist = filelist, which.lines = 1000)
#'
#' # option 2: specifiying a directory
#' fcs2 <- read.FCSset(directory = "../data", pattern = ".LMD", num.threads = 4)
#' }

read.FCSset <- function(filelist = NULL, directory = getwd(), pattern = ".fcs$", 
                        events = "all", dataset = 1, num.threads = NULL){
  set.seed(333)
  if(is.null(filelist)) print(filelist <- list.files(path = directory, pattern = pattern, full.names = T))
  
  if(Sys.info()[1] == "Windows") require(parallelsugar) #devtools::github_install('nathanvan/parallelsugar')
  if(events == "all"){
  	fcs <- as.flowSet(read.ncdfFlowSet(filelist, emptyValue = FALSE, transformation = FALSE, truncate_max_range = FALSE, 
                                which.lines = NULL, dataset = dataset, mc.cores = num.threads))
  	}else{
	fcs <- as.flowSet(read.ncdfFlowSet(filelist, emptyValue = FALSE, transformation = FALSE, truncate_max_range = FALSE, 
                                which.lines = events, dataset = dataset, mc.cores = num.threads))
  	}
  
  return(fcs)
}
