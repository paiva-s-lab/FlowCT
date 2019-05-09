#' reduce.FCS
#'
#' This function reduce the number of events for a single FCS file or multiple FCS files contained in a given folder.
#' @param file_or_directory FCS filename or directory where to search for FCS files. Default = "."
#' @param keep_n_events Number of events to keep in each FCS files. Default = 100000
#' @param output_folder Folder name for storing new reduced FCS files. Default = "reduced"
#' @param output_suffix Suffix to be added to reduced FCS files. Default = "_red"
#' @keywords reduce
#' @keywords events
#' @export
#' @examples
#' \dontrun{reduce.FCS()}

reduce.FCS <- function(file_or_directory = ".", keep_n_events = 100000, output_folder = "reduced",
                       output_suffix = "_red"){

  if(fs::is_file(file_or_directory)){
    filelist <- file_or_directory
  }else{
    cat("FCS files within the selected folder:\n")
    print(filelist <- list.files(path = file_or_directory, pattern = ".fcs$", full.names = F))
    cat("\n")
  }

  dir.create(output_folder)

  for(file in filelist){
    cat(paste0("Processing: ", file, "\n"))

    aux <- flowCore::read.FCS(file)
    if(dim(aux@exprs)[1] < keep_n_events){
      cat(paste0(file, " has a lower number of events that keep_n_events so it'll be only copied.\n"))
      flowCore::write.FCS(aux, paste0(output_folder, "/", gsub(".fcs","", file), output_suffix, ".fcs"))
    }else{
      tmp_exprs <- aux@exprs[sample(nrow(aux@exprs), keep_n_events),] #reduce events number

      aux@parameters@data$desc <- ifelse(is.na(aux@parameters@data$desc), aux@parameters@data$name, aux@parameters@data$desc) #replace NA values with CD names
      aux@parameters@data$desc <- sapply(aux@parameters@data$desc, function(x) strsplit(x, ":")[[1]][1]) #separe maker:antibody (eg. CCR4:PE-Cy7-A)

      #generate a new FCS
      aux <- methods::new("flowFrame", exprs  = tmp_exprs, parameters = aux@parameters, description = aux@description)
      flowCore::write.FCS(aux, paste0(output_folder, "/", gsub(".fcs","", file), output_suffix, ".fcs"))
    }
  }
  rm(aux, tmp_exprs)
}

remove.doublets <- function(flowFrame, d1="FSC-A", d2="FSC-H", w=NULL,silent=TRUE){
  #calculate the ratios
  ratio <- flowCore::exprs(flowFrame)[,d1] / (1+ flowCore::exprs(flowFrame)[,d2])

  #define the region that is accepted
  r <- stats::median(ratio)
  if(is.null(w)){ w <- 2*stats::sd(ratio) }
  if(!silent){
    print(r)
    print(w)
  }

  #make selection
  selection <- which(ratio < r+w)
  return(flowFrame[selection,])
}


#' qc.and.removeDoublets
#'
#' This function performs a double-step procedure in all FCS files cointained within a folder:\enumerate{
#'    \item It executes a quality control based in \href{http://bioconductor.org/packages/release/bioc/vignettes/flowAI/inst/doc/flowAI.html}{\code{\link[flowAI:flow_auto_qc]{flowAI::flow_auto_qc()}}},
#'    \item Remove doublets according criteria described \href{https://github.com/LennonLab/flow-cytometry/blob/fcf09fc7b4943a386864de9b8ee43ce4e1c8d34d/bin/support_functions.R}{here}}.
#' @param directory Directory with FCS files to process. Default = "."
#' @param reduction_computed Logical whether reduction previously computed with \code{\link[FlowCT:reduce.FCS]{FlowCT::reduce.FCS()}}. Default = \code{TRUE}
#' @param reduction_suffix Suffix to consider in reduced files (normaly it would be "_red" if it is coming from \code{\link[FlowCT:reduce.FCS]{FlowCT::reduce.FCS()}}). Default = "_red"
#' @param output_suffix Suffix to be added to final FCS files. Default = "_preprocessed"
#' @param output_folder Folder name for storing final FCS files. Default = "results_HQsinglets"
#' @keywords quality control
#' @keywords doublets
#' @keywords remove
#' @export
#' @examples
#' \dontrun{qc.and.removeDoublets()}

qc.and.removeDoublets <- function(directory = ".",  reduction_computed = TRUE, reduction_suffix = "_red",
                                  output_suffix = "_preprocessed", output_folder = "results_HQsinglets"){
  setwd(directory)
  dir.create(output_folder) #create the output folder

  tryCatch( #handle errors
    {if(reduction_computed){
      filelist <- list.files(pattern = paste0(reduction_suffix, ".fcs$"))
    }else{
      filelist <- list.files(pattern = ".fcs$")
    }

      suppressWarnings(flowAI::flow_auto_qc(filelist, ChExcludeFS = c("FSC-H", "SSC-H","FSC-A", "SSC-A"), fcs_highQ = "_HQ",
                                    fcs_QC = FALSE, output = 0, html_report = F, folder_results = F, mini_report = "miniQC"))

      ## Doublet removal
      if(reduction_computed){
        filelist <- list.files(pattern = paste0(reduction_suffix, "_HQ.fcs$"))
      }else{
        filelist <- list.files(pattern = "_HQ.fcs$")
      }


      for(file in filelist){
        flowFrame <- flowCore::read.FCS(paste(file,sep="/"))
        flowFrame_d <- remove.doublets(flowFrame, d1 = "FSC-A", d2 = "FSC-H")

        if(reduction_computed){
          flowCore::write.FCS(flowFrame_d, paste0(output_folder, "/", gsub(paste0(reduction_suffix, "_HQ.fcs$"), "", file), output_suffix, ".fcs"))
        }else{
          flowCore::write.FCS(flowFrame_d, paste0(output_folder, "/", gsub("_HQ.fcs", "", file), output_suffix, ".fcs"))
        }
      }

      ## calculate ratios for corrected samples
      if(reduction_computed){
        filelist <- list.files(pattern = paste0(reduction_suffix, ".fcs$"))
      }else{
        filelist <- list.files(pattern = ".fcs$")
      }

      for(i in filelist){
        if(reduction_computed){
          n_aux <- gsub("red.fcs", "", i)
        }else{
          n_aux <- gsub(paste0(reduction_suffix, "_red_HQ.fcs$"), "", i)
        }

        r_aux <- flowCore::read.FCS(i)
        red <- as.numeric(dim(r_aux@exprs)[1])

        f_aux <- flowCore::read.FCS(paste0(output_folder, "/", list.files(path = output_folder))[grep(n_aux, list.files(path = output_folder), fixed = T)])
        final <- as.numeric(dim(f_aux@exprs)[1])

        if(final/red < 0.7){
          print(paste0("WARNING! > ", i, " has lost some much cells (more that 30%) in the QC and doublets removal steps, consider to review it!"))
        }
      }

      qct <- utils::read.table("miniQC.txt", header = T, sep = "\t", row.names = 1,
                        colClasses = c(rep("character", 2), rep("NULL", 4)))
      gridExtra::grid.arrange(gridExtra::tableGrob(qct, cols = c("# initial events", "% deleted events")))

      invisible(file.remove(list.files(pattern = "_HQ.fcs$|miniQC")))},

    error = function(e) {invisible(file.remove(list.files(pattern = "_HQ.fcs$|miniQC"))); print("An ERROR has occurred!")})
}
