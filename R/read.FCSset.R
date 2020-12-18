#' Load multiple FCS files
#'
#' It reads multiple FCS files (in a computationaly-distributed way) contained in a specific folder or given through a vector. Important note: all files should have an identicar header (i.e., same name and markers order), if not only those common will be read (you can use \code{\link[FlowCT:unify.FCSheaders]{FlowCT::unify.FCSheaders()}} for doing that).
#' @param filelist A vector with full path of FCS files to be read, commonly generated through \code{\link[base:list.files]{base::list.files()}}. If \code{NULL}, this file list will be generated as indicated below.
#' @param directory If \code{filelist = NULL}, those files stored in this location will be read. Default = \code{getwd()} (current directory).
#' @param pattern Pattern for reading files within \code{directory}. Default = \code{"fcs"}.
#' @param events Numeric vector indicating how many events to read in each FCS file (same behaviour as \code{\link[flowCore:read.FCS]{flowCcre::read.FCS()}}). Default = \code{NULL}, i.e., all events will be read.
#' @param dataset Specify data segment to be read in the FCS file (same behaviour as \code{\link[flowCore:read.FCS]{flowCcre::read.FCS()}}). Default = \code{1}.
#' @param num.threads Number of threads for reading the FCS files. For Windows users, you should to install \code{parallelsugar} package from (\href{https:https://github.com/nathanvan/parallelsugar}{Nathanvan's GitHub}). Defult = \code{NULL} (i.e., not parallelization).
#' @keywords FCS reading
#' @keywords FCS parallel
#' @export
#' @importFrom ncdfFlow as.flowSet read.ncdfFlowSet
#' @importFrom flowCore fsApply exprs sampleNames
#' @importFrom stats complete.cases
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
  if(is.null(filelist)) filelist <- list.files(path = directory, pattern = pattern, full.names = T)
  print(basename(filelist))

  if(Sys.info()[1] == "Windows") require(parallelsugar) #devtools::github_install('nathanvan/parallelsugar')
  suppressMessages(fcs <- as.flowSet(read.ncdfFlowSet(filelist, emptyValue = FALSE, transformation = FALSE, truncate_max_range = FALSE,
                                                      which.lines = NULL, dataset = dataset, mc.cores = num.threads)))

  ## downsample FCS
  if(events != "all"){
    for(i in fcs@phenoData@data$name){
      if(nrow(fcs[[i]]) < events){
        warning("Filename: ", i, " has a lower number of events than specified through 'events', all will be read!", call. = F, immediate. = T)
      }else{
        fcs[[i]] <- fcs[[i]][sample(nrow(fcs[[i]]), events)]
      }
    }
  }

  ## check NA events because previous flow cytometer exporting errors
  wrong <- fsApply(fcs, function(x) sum(apply(exprs(x), 2, is.na)))[,1]
  wrong <- names(wrong)[wrong != 0]

  for(i in wrong){
    aux <- fcs[[i]]
    if(sum(apply(exprs(aux), 2, is.na))!= 0)
      cat("File", i, ">", sum(apply(exprs(aux), 2, is.na))/ncol(aux), "NA events deleted (flow cytometer exporting errors)\n")
  exprs(aux) <- exprs(aux)[complete.cases(exprs(aux)),]
  fcs[[i]] <- aux
  }

  return(fcs)
}
