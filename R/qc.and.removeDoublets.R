#' qc.and.removeDoublets
#'
#' This function performs a quality control based on a double-step procedure in all FCS files contained within a \code{FCS.SE} object or a folder specified by user:\enumerate{
#'    \item It executes a quality control based in \href{http://bioconductor.org/packages/release/bioc/vignettes/flowAI/inst/doc/flowAI.html}{\code{flowAI::flow_auto_qc()}}},
#'    \item Remove doublets according criteria described \href{https://github.com/LennonLab/flow-cytometry/blob/fcf09fc7b4943a386864de9b8ee43ce4e1c8d34d/bin/support_functions.R}{here}}.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}. If \code{NULL} (default), this function will search for FCS files within current working directory.
#' @param filelist A vector with full path of FCS files to be read, commonly generated through \code{\link[base:list.files]{base::list.files()}}. If \code{NULL}, this file list will be generated as indicated below.
#' @param directory If \code{filelist = NULL}, those files stored in this location will be read. Default = \code{getwd()} (current directory).
#' @param pattern Pattern for reading files within \code{directory}. Default = \code{"fcs"}.
#' @param physical.markers Vector with physical markers (i.e., FCS-A, SSC-H, etc.). Mandatory for perfoming quality control.
#' @param output.folder Folder name for storing final FCS files. Default = \code{NULL} (current directory).
#' @param output.suffix Suffix added to new generated FCS files (in case not being working with a \code{FCS.SE} object). Default = \code{qc}.
#' @param return.idx Logical indicating whether output is a new \code{FCS.SE} object without low-quality events or an index with this positions (without removing them from original \code{FCS.SE}). This option is only available if the input is a \code{FCS.SE} object. Default = \code{F} (\code{FCS.SE} output).
#' @keywords quality control
#' @keywords doublets removal
#' @keywords remove low quality events
#' @export
#' @importFrom flowCore exprs
#' @return Final output is a new \code{FCS.SE} object (if FCS.SE input) without these events that do not pass the quality control or new FCS files with the suffix \code{qc} with these low quality events removed. If user specifies \code{return.idx = T}, output would be a vector with all low-quality events positions.
#' @return A table with percentaje of removed events for each FCS will be also shown in the terminal (those files with more than 30% of events removed will have a '!' signal).
#' @examples
#' \dontrun{
#' # option 1: output a new FCS.SE object with low-quality events removed.
#' fcs_se_qc <- qc.and.removeDoublets(fcs.SE = fcs_se, physical.markers = c("FSC_A", "FSC_H", "SSC_A", "SSC_H"))
#'
#' # option 2: output an index with low quality events
#' idx_qc <- qc.and.removeDoublets(fcs.SE = fcs_se, physical.markers = c("FSC_A", "FSC_H", "SSC_A", "SSC_H"), return.idx = T)
#' 
#' # option 3: working with all FCS files within a folder
#' idx_qc <- qc.and.removeDoublets(directory = "../data/", physical.markers = c("FSC_A", "FSC_H", "SSC_A", "SSC_H"), output.folder = "HQ_files")
#' }

qc.and.removeDoublets <- function(fcs.SE = NULL, filelist = NULL, directory = getwd(), pattern = "fcs", physical.markers, output.folder = NULL, output.suffix = "qc", return.idx = F){
  if(!is.null(fcs.SE)){
    fcs <- as.flowSet.SE(fcs.SE, assay.i = "raw")
    filenames <- fcs@phenoData@data$name
    
    idx <- losses <- c()
    for(file in filenames){
      fcs[[file]]@description$FILENAME <- file #requirements for flow_auto_qc
      
      invisible(utils::capture.output(idx1 <- flowAI::flow_auto_qc(fcs[[file]], ChExcludeFS = physical.markers, fcs_QC = FALSE, output = 3, html_report = F, folder_results = F, mini_report = F)[[1]]))
      
      # doublet removal
      FSCA <- grep("FS.*.A", physical.markers, value = T)
      FSCH <- grep("FS.*.H", physical.markers, value = T)
      
      ratio <- exprs(fcs[[file]])[,FSCA]/(1+ exprs(fcs[[file]])[,FSCH]) #calculate the ratios
      r <- stats::median(ratio)
      w <- 2*stats::sd(ratio) 
      idx2 <- as.vector(which(ratio > r+w)) #define the region that will be removed
      
      idx <- c(idx, rownames(exprs(fcs[[file]]))[c(idx1, idx2)])
      losses <- append(losses, round(length(c(idx1, idx2))/nrow(fcs[[file]])*100, digits = 2))
      names(losses)[length(losses)] <- file
    }
    print(knitr::kable(as.data.frame(losses), col.names = "removed events (%)"))
    if(return.idx) return(idx) else return(fcs.SE[,-match(idx, colnames(fcs.SE))])
  }else{
    if(!is.null(output.folder)) dir.create(output.folder) else output.folder <- directory 
    if(is.null(filelist)) filelist <- list.files(path = directory, pattern = pattern, full.names = T)
    
    for(file in filelist){
      print(file)
      invisible(utils::capture.output(file_q <- flowAI::flow_auto_qc(file, ChExcludeFS = physical.markers, fcs_QC = FALSE, output = 1, html_report = F, folder_results = F, mini_report = "miniQC")))
      
      # doublet removal
      FSCA <- grep("FS.*.A", physical.markers, value = T)
      FSCH <- grep("FS.*.H", physical.markers, value = T)
      
      ratio <- exprs(file_q)[,FSCA]/(1+ exprs(file_q)[,FSCH]) #calculate the ratios
      r <- stats::median(ratio)
      w <- 2*stats::sd(ratio) 
      file_d <- file_q[which(ratio < r+w),] #define the region that is accepted
      
      extension <- strsplit(basename(file), "\\.")[[1]][2]
      filename <- strsplit(basename(file), "\\.")[[1]][1]
      flowCore::write.FCS(file_d, paste0(output.folder, "/", filename, ".", output.suffix, ".", extension))
    }
    # notify removed evens
    qct <- utils::read.table("miniQC.txt", header = T, sep = "\t", colClasses = c("character", rep("numeric", 2), rep("NULL", 4)))
    qct$warning <- ifelse(qct$X..anomalies >= 30, "(!)", "")
    print(knitr::kable(qct, col.names = c("filenames", "# initial events", "% deleted events", "Warning!")))
  }
}
