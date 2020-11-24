#' Perform quality control and remove doublets
#'
#' This function performs a quality control based on a double-step procedure in all FCS files contained within a \code{fcs.SCE} object or a folder specified by user:\enumerate{
#'    \item It executes a quality control based in \href{http://bioconductor.org/packages/release/bioc/vignettes/flowAI/inst/doc/flowAI.html}{\code{flowAI::flow_auto_qc()}},
#'    \item Remove doublets according criteria described \href{https://github.com/LennonLab/flow-cytometry/blob/fcf09fc7b4943a386864de9b8ee43ce4e1c8d34d/bin/support_functions.R}{here}.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT.v2:fcs.SCE]{FlowCT.v2::fcs.SCE()}}. If \code{NULL} (default), this function will search for FCS files within current working directory.
#' @param filelist A vector with full path of FCS files to be read, commonly generated through \code{\link[base:list.files]{base::list.files()}}. If \code{NULL}, this file list will be generated as indicated below.
#' @param directory If \code{filelist = NULL}, those files stored in this location will be read. Default = \code{getwd()} (current directory).
#' @param pattern Pattern for reading files within \code{directory}. Default = \code{"fcs"}.
#' @param physical.markers Vector with physical markers (i.e., FCS-A, SSC-H, etc.). Mandatory for perfoming quality control.
#' @param output.folder Folder name for storing final FCS files. Default, the same of \code{directory}.
#' @param output.suffix Suffix added to new generated FCS files (in case not being working with a \code{fcs.SCE} object). Default = \code{qc}.
#' @param return.idx Logical indicating whether output is a new \code{fcs.SCE} object without low-quality events or an index with this positions (without removing them from original \code{fcs.SCE}). This option is only available if the input is a \code{fcs.SCE} object. Default = \code{FALSE} (\code{fcs.SCE} output).
#' @param return.fcs Logical indicating whether new FCS file should be generated without low quality events. They will stored in specified \code{output.folder} with the corresponding \code{output.suffix}.
#' @param keep.QCfile Final QC text report should be kept? Attention, this file does not contain removed doublets but only those events removed through {\code{flowAI::flow_auto_qc()}}, all events are showed trough terminal output. Default = \code{FALSE}.
#' @keywords quality control
#' @keywords doublets removal
#' @keywords remove low quality events
#' @export
#' @importFrom flowCore exprs write.FCS
#' @importFrom premessa as_flowFrame
#' @importFrom flowAI flow_auto_qc
#' @importFrom stats median sd
#' @importFrom utils capture.output read.table
#' @importFrom knitr kable
#' @return Final output is a new \code{fcs.SCE} object (if \code{fcs.SCE} input) without these events that do not pass the quality control or new FCS files with the suffix \code{"qc"} with these low quality events removed. If user specifies \code{return.idx = TRUE}, output would be a vector with all low-quality events positions.
#' @return A table with percentaje of removed events for each FCS will be also shown in the terminal (those files with more than 30% of events removed will have a '!' signal).
#' @examples
#' \dontrun{
#' # option 1: output a new fcs.SCE object with low-quality events removed.
#' fcs.SCE_qc <- qc.and.removeDoublets(fcs.SCE = fcs, 
#'      physical.markers = c("FSC_A", "FSC_H", "SSC_A", "SSC_H"))
#'
#' # option 2: output an index with low quality events
#' idx_qc <- qc.and.removeDoublets(fcs.SCE = fcs, 
#'      physical.markers = c("FSC_A", "FSC_H", "SSC_A", "SSC_H"), return.idx = T)
#' 
#' # option 3: working with all FCS files within a folder
#' idx_qc <- qc.and.removeDoublets(directory = "../data/", output.folder = "HQ_files",
#'      physical.markers = c("FSC_A", "FSC_H", "SSC_A", "SSC_H"))
#' }

qc.and.removeDoublets <- function(fcs.SCE = NULL, filelist = NULL, directory = getwd(), pattern = "fcs", physical.markers, output.folder = directory, output.suffix = "qc", return.idx = FALSE, return.fcs = FALSE, keep.QCfile = F){
  if(!is.null(fcs.SCE)){
    fcs <- as.flowSet.SE(fcs.SCE, assay.i = "raw")
    filenames <- fcs@phenoData@data$name
    
    idx <- losses <- c()
    for(file in filenames){
      fcs[[file]]@description$FILENAME <- file #requirements for flow_auto_qc
      
      invisible(capture.output(idx1 <- flow_auto_qc(fcs[[file]], ChExcludeFS = physical.markers, fcs_QC = FALSE, output = 3, html_report = F, folder_results = F, mini_report = "miniQC")[[1]]))
      
      # doublet removal
      FSCA <- grep("FS.*.A", physical.markers, value = T)
      FSCH <- grep("FS.*.H", physical.markers, value = T)
      
      exprstmp <- exprs(fcs[[file]])

      ratio <- exprstmp[,grep(FSCA, colnames(exprstmp))]/(1+ exprstmp[,grep(FSCH, colnames(exprstmp))]) #calculate the ratios
      r <- median(ratio)
      w <- 2*sd(ratio) 
      idx2 <- as.vector(which(ratio > r+w)) #define the region that will be removed
      
      idx <- append(idx, rownames(exprstmp)[c(idx1, idx2)])
      losses <- append(losses, round(length(c(idx1, idx2))/nrow(fcs[[file]])*100, digits = 2))
      names(losses)[length(losses)] <- file

      # generate new high QC FCS
      if(return.fcs){
        extension <- tolower(tail(strsplit(filenames[1], split="\\.")[[1]], n = 1))
        idx_fcs <- unique(rownames(exprstmp)[c(idx1, idx2)])

        exp_matrix <- exprstmp[-match(idx_fcs, rownames(exprstmp)),]
        colnames(exp_matrix) <- sapply(fcs.SCE@metadata$raw_channel_names, function(x) strsplit(x, split = ":")[[1]][1])
        file_return <- as_flowFrame(exp_matrix)
        # file_return <- flowCore::flowFrame(exp_matrix)
        file_return@parameters@data$desc <- sapply(fcs.SCE@metadata$raw_channel_names, function(x) strsplit(x, split = ":")[[1]][2])
                                                   
        write.FCS(file_return, paste0(output.folder, "/", gsub(extension, "", file), output.suffix, ".", extension))
      }
    }

    if(return.fcs) cat("New generated FCSs without low quality events are stored in:", output.folder)
    print(kable(as.data.frame(losses), col.names = "removed events (%)"))
    fcs.SCE@metadata$lowQC_events <- idx

    if(!keep.QCfile) invisible(file.remove(list.files(pattern = "miniQC")))

    if(return.idx) return(idx) else return(fcs.SCE[,-match(idx, colnames(fcs.SCE))])
  }else{
    if(is.null(filelist)) filelist <- list.files(path = directory, pattern = pattern, full.names = T)
    
    for(file in filelist){
      print(file)
      invisible(capture.output(file_q <- flow_auto_qc(file, ChExcludeFS = physical.markers, fcs_QC = FALSE, output = 1, html_report = F, folder_results = F, mini_report = "miniQC")))
      
      # doublet removal
      FSCA <- grep("FS.*.A", physical.markers, value = T)
      FSCH <- grep("FS.*.H", physical.markers, value = T)
      
      exprstmp <- exprs(fcs[[file]])

      ratio <- exprstmp[,grep(FSCA, colnames(exprstmp))]/(1+ exprstmp[,grep(FSCH, colnames(exprstmp))]) #calculate the ratios
      r <- median(ratio)
      w <- 2*sd(ratio) 
      file_d <- file_q[which(ratio < r+w),] #define the region that is accepted
      
      extension <- tolower(tail(strsplit(basename(file), "\\.")[[1]], n = 1))
      filename <- strsplit(basename(file), "\\.")[[1]][1]
      write.FCS(file_d, paste0(output.folder, "/", filename, ".", output.suffix, ".", extension))
    }
    # notify removed evens
    qct <- read.table("miniQC.txt", header = T, sep = "\t", colClasses = c("character", rep("numeric", 2), rep("NULL", 4)))
    qct$warning <- ifelse(qct$X..anomalies >= 30, "(!)", "")
    print(kable(qct, col.names = c("filenames", "# initial events", "% deleted events", "Warning!")))

    if(!keep.QCfile) invisible(file.remove(list.files(pattern = "miniQC")))
  }
}
