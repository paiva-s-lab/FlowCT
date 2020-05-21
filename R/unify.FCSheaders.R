#' unify.FCSheaders 
#'
#' It checks if FCS files within a specific folder, or indicated in a vector, have the same header (i.e., same channels nomenclature and order) and offer the possibility to unify them creating new ones.
#' @param filelist A vector with full path of FCS files to be read, commonly generated through \code{\link[base:list.files]{base::list.files()}}. If \code{NULL}, this file list will be generated as indicated below.
#' @param directory If \code{filelist = NULL}, those files stored in this location will be read. Default = \code{getwd()} (current directory).
#' @param pattern Pattern for reading files within \code{directory}. Default = \code{"fcs"}.
#' @param dataset Specify data segment to be read in the FCS file (same behaviour as \code{\link[flowCore:read.FCS]{flowCore::read.FCS()}}). Default = \code{1}.
#' @param fix Logical indicating if user wants to perform the unification of FCS's headers or only display header frecuencies. If \code{TRUE} (default value), you should to specify the chosen frequency.
#' @param select.freq Of diplayed frequencies (row-numered), indicate which is your option for unifiying FCS files. Default = \code{1} (the more frequent). whether
#' @param events Numeric vector indicating how many events to read in each FCS file (same behaviour as \code{\link[flowCore:read.FCS]{flowCore::read.FCS()}}), those are the number of events for new and unified generated FCS files. Default = \code{NULL}, i.e., all events will be read.
#' @param view Logical indicating whether table with frequencies should be displayed in the terminal or in a external table. Default = \code{FALSE}.
#' @keywords unify FCS headers
#' @keywords differing FCS files
#' @return The final output if \code{fix = F} is a table with three columns: \code{names} for channels, \code{freq} for frequency of apparition of these channel names and \code{length} for number of channels within that frequency.
#' @return If \code{fix = T}, those FCS files with a distinct header from selected frequency will be renamed/reordered, added with the suffix \code{fixed.fcs} and the original ones (unchanged) will be stored in a new \code{original_files} folder. In the case a FCS would have a different number of channels, it will moved to a new folder called \code{discarded_files_because_diffs} and dicarded from downstream analysis (this will be changed in the future).
#' @export
#' @examples
#' \dontrun{
#' ## detect header's frequencies
#' unify.FCSheaders(directory = "../data/", pattern = "fcs", fix = F)
#'  #                                                                            names freq length
#'  # 1 FSC-H, SSC-H, FITC-A, PE-A, PerCP-Cy5-5-A, APC-A, APC-H7-A, V500-A, V450-A       12     10
#'  # 2 FSC-H, SSC-H, FITC-A, PE-A, PerCP-C5-A, APC-A, APC-H7-A, Pac-Orange, Pac-Blue    10     10
#' ## unify FCS's headers
#' unify.FCSheaders(directory = "../data/", pattern = "fcs", fix = T, select.freq = 2)
#' }

unify.FCSheaders <- function(filelist = NULL, directory = getwd(), pattern = ".fcs$", dataset = 1, fix = T, select.freq = 1, events = "all", view = F){
  
  if(is.null(filelist)) filelist <- list.files(path = directory, pattern = pattern, full.names = T)
  if(!is.null(filelist)) directory <- ""
  if(events == "all") events <- NULL

  
  ## test frequency for each combination of channel/markers
  mnames <- c()
  for(i in 1:length(filelist)){
    aux <- flowCore::read.FCS(filelist[i], dataset = dataset, emptyValue = FALSE, transformation = FALSE, truncate_max_range = FALSE, which.lines = 1)
    mnames <- append(mnames, paste(colnames(aux), collapse = ", "))
  }
  
  tm <- data.table::melt(table(mnames))
  colnames(tm) <- c("names", "freq")
  tm$length <- sapply(tm$names, function(x) length(strsplit(as.character(x), ", ")[[1]]))
  tm <- tm[order(tm$freq, decreasing = T),]
  rownames(tm) <- NULL
  
  if(fix){
    ## choose frequency
    if(nrow(tm) == 1){
      return(cat("All files are correctly and uniformly named!\n"))
    }else{
      cat("--------------------------\nTable with frequency of differently named FCSs:\n")
      if(view) utils::View(tm) else print(tm)
      
      if(select.freq == "ask"){
        select.freq <- as.numeric(readline("Select a number for establishing a common marker order: "))
        o <- strsplit(as.character(tm$names[select.freq]), ", ")[[1]]
      }else{
        o <- strsplit(as.character(tm$names[select.freq]), ", ")[[1]]
      }
    } 
    
    ## detect files that not keep more abundant order and reorder/rename them
    pb <- progress::progress_bar$new(total = length(filelist), format = "\nCorrecting divergent files [:bar] :percent eta: :eta")
    
    for(i in 1:length(filelist)){
      pb$tick()
      # diff_files <- c()
      aux <- flowCore::read.FCS(filelist[i], dataset = dataset, emptyValue = FALSE, transformation = FALSE, truncate_max_range = FALSE, which.lines = events)
      
      if(!identical(o, colnames(aux))){
        if(length(colnames(aux)) != length(o)){
          cat("\nFile", filelist[i], "has a some different markers or channels, it cannot be included in the analysis. Moved to 'discarded_files_because_diffs/'\n")
          dir.create(paste0(directory, "discarded_files_because_diffs"))
          suppressMessages(filesstrings::file.move(filelist[i], paste0(directory, "discarded_files_because_diffs")))
          # diff_files <- append(diff_files, filelist[i])
        }else{
          dir.create(paste0(directory, "original_files"), showWarnings = F) #move files
          suppressMessages(filesstrings::file.move(filelist[i], paste0(directory, "original_files")))
          
          aux@parameters@data <- aux@parameters@data[match(o, aux@parameters@data$name),] #reorder
          colnames(aux) <- o
          
          extension <- strsplit(filelist[i], "\\.")
          extension <- extension[[1]][length(extension[[1]])]
          flowCore::write.FCS(aux, filename = paste0(gsub(extension, "", filelist[i]), "fixed.", extension))
        }      
      }
      Sys.sleep(1/100)
    }
  }else{
    if(view) utils::View(tm) else print(tm)
  }
}