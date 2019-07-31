#' reduce.FCS
#'
#' This function reduce the number of events for a single FCS file or multiple FCS files contained in a given folder.
#' @param file_or_directory FCS filename or directory where to search for FCS files. Default = "." (current directory)
#' @param keep_n_events Number of events to keep in each FCS files. Default = 100000
#' @param output_suffix Suffix to be added to reduced FCS files. Default = "_red"
#' @param output_folder Folder name for storing new reduced FCS files. Default = "reduced"
#' @keywords reduce
#' @keywords events
#' @export
#' @importFrom fs is_file
#' @examples
#' \dontrun{reduce.FCS()}

reduce.FCS <- function(file_or_directory = ".", keep_n_events = 100000, output_suffix = "_red",
                       output_folder = "reduced"){
  if(is_file(file_or_directory)){
    filelist <- file_or_directory
  }else{
    cat("FCS files within the selected folder:\n")
    print(filelist <- list.files(path = file_or_directory, pattern = ".fcs$", full.names = F))
    cat("\n")
  }
  
  dir.create(output_folder)
  
  for(file in filelist){
    cat(paste0("Processing: ", file, "\n"))
    
    aux <- read.FCS(file)
    if(dim(aux@exprs)[1] < keep_n_events){
      cat(paste0(file, " has a lower number of events that keep_n_events so it'll be only copied.\n"))
      write.FCS(aux, paste0(output_folder, "/", gsub(".fcs","", file), output_suffix, ".fcs"))
    }else{
      tmp_exprs <- aux@exprs[sample(nrow(aux@exprs), keep_n_events),] #reduce events number
      
      aux@parameters@data$desc <- ifelse(is.na(aux@parameters@data$desc), aux@parameters@data$name, aux@parameters@data$desc) #replace NA values with CD names
      aux@parameters@data$desc <- sapply(aux@parameters@data$desc, function(x) strsplit(x, ":")[[1]][1]) #separe maker:antibody (eg. CCR4:PE-Cy7-A)
      
      #generate a new FCS
      aux <- new("flowFrame", exprs  = tmp_exprs, parameters = aux@parameters, description = aux@description)
      write.FCS(aux, paste0(output_folder, "/", gsub(".fcs","", file), output_suffix, ".fcs"))
    }
  }
  rm(aux, tmp_exprs)
  setwd(output_folder)
}

                                         
#' qc.and.removeDoublets
#'
#' This function performs a double-step procedure in all FCS files cointained within a folder:\enumerate{
#'    \item It executes a quality control based in \href{http://bioconductor.org/packages/release/bioc/vignettes/flowAI/inst/doc/flowAI.html}{\code{flowAI::flow_auto_qc()}}},
#'    \item Remove doublets according criteria described \href{https://github.com/LennonLab/flow-cytometry/blob/fcf09fc7b4943a386864de9b8ee43ce4e1c8d34d/bin/support_functions.R}{here}}.
#' @param directory Directory with FCS files to process. Default = "." (current directory)
#' @param reduction_suffix Suffix added previously to reduced files trough \emph{Rscript consolidate_and_reduction.rscript} or \code{\link[FlowCT:reduce.FCS]{FlowCT::reduce.FCS()}}. If \code{NULL}, you're indicating previous reduction has not been performed. Default = ".red"
#' @param output_suffix Suffix to be added to final FCS files. Default = ".preprocessed"
#' @param output_folder Folder name for storing final FCS files. Default = "results_preprocessing"
#' @keywords quality control
#' @keywords doublets
#' @keywords remove
#' @export
#' @importFrom flowAI flow_auto_qc
#' @importFrom ggplot2 ggsave
#' @importFrom utils read.table
#' @import flowCore
#' @import gridExtra
#' @examples
#' \dontrun{qc.and.removeDoublets()}

qc.and.removeDoublets <- function(directory = ".", reduction_suffix = ".red",
                                  output_suffix = ".preprocessed", output_folder = "results_preprocessing"){
  
  lapply(c("flowAI", "flowCore", "gridExtra", "ggplot2"), require, character.only = TRUE)
  
  setwd(directory)
  dir.create(output_folder) #create the output folder
  
  if(!is.null(reduction_suffix)){
    filelist <- list.files(pattern = paste0(reduction_suffix, ".fcs$"))
  }else{
    filelist <- list.files(pattern = ".fcs$")
  }
    
  for(file in filelist){
    cat(paste0("Processing: ", file, "\n")) 
    sink("tmp") #to turn off stdout (other options don't work)
    file_q <- suppressWarnings(flow_auto_qc(file, ChExcludeFS = c("FSC-H", "SSC-H","FSC-A", "SSC-A"),
                                  fcs_QC = FALSE, output = 1, html_report = F, folder_results = F, mini_report = "miniQC"))
    sink();invisible(file.remove(list.files(pattern = "tmp")))
    
    ## doublet removal
    ratio <- exprs(file_q)[,"FSC-A"]/(1+ exprs(file_q)[,"FSC-H"]) #calculate the ratios
    r <- median(ratio)
    w <- 2*sd(ratio) 
    file_d <- file_q[which(ratio < r+w),] #define the region that is accepted

    ## calculate ratios for corrected samples
    orig <- read.FCS(file)
    orig <- as.numeric(dim(orig@exprs)[1])
    red <- as.numeric(dim(file_d@exprs)[1])
      
    if(red/orig < 0.7){
      cat("WARNING! >", i, "has lost some much cells (more that 30%) in the QC and doublets removal steps, consider to review it!", )
    }
  }
    
  cat("------------------------------\nFinal QC and remove doublets result:")  
  qct <- read.table("miniQC.txt", header = T, sep = "\t", 
                    colClasses = c("character", rep("numeric", 2), rep("NULL", 4)))
  qct$warning <- ifelse(qct$X..anomalies >= 30, "(!)", "")
  print(knitr::kable(qct, col.names = c("filenames", "# initial events", "% deleted events", "Warning!")))
  cat("\n")

  write.FCS(file_d, paste0(output_folder, "/", gsub(".fcs$", "", file), output_suffix, ".fcs"))
  setwd(output_folder)
  invisible(file.remove(list.files(pattern = "miniQC")))
}
