#' sub.samples
#'
#' It generates reduce randomly the number of events of a FCS.SE object (it computes this reduction for each FCS file separatelly inside this object). It can generate a new reduced FCS.SE object or a simple index position for removing events.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param subsampling Number of events to keep. Default = 1000.
#' @param index Logical indicating if returning an FCS.SE object or a index vector. Default = \code{FALSE}.
#' @keywords subsampling
#' @export
#' @examples
#' \dontrun{
#' # option 1
#' fcs_se_red <- sub.samples(data = fcs_se, subsampling = 1000)
#'
#' # option 2
#' keep_red <- sub.samples(data = fcs_se, subsampling = 1000, index = T)
#' fcs_se_red <- fcs_se[,keep_red]
#' }

sub.samples <- function (fcs.SE, subsampling = 1000, index = F){
  set.seed(333)
  pb <- progress::progress_bar$new(total = length(unique(fcs.SE$filename)), format = "Random subsampling [:bar]")
  
  sub_idx <- vector()
  for(i in unique(fcs.SE$filename)) {
    pb$tick()
    
    aux <- colnames(fcs.SE[,fcs.SE$filename == i])
    if(length(aux) < subsampling){
      cat("Attention!: filename", i, "has a lower number of events, reduction won't be computed.\n")
      sub_idx <- append(sub_idx, aux)
    }else{
      sub_idx <- append(sub_idx, aux[sample(length(aux), subsampling)])    
    }
    
    Sys.sleep(1/10)
  }
  if(index) return(sub_idx) else return(fcs.SE[,sub_idx])
}

                                                                                        
#' export.metaFCS
#'
#' It creates a FCS file containing all analysis incorpored to colData(fcs.SE) as well as the dimensional reduction coordinates calculated.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param dr.object Object created with \code{\link[FlowCT:dim.reduction]{FlowCT::dim.reduction()}} function or a table combining DR and experimental metadata information.
#' @param output.name Name for generated FCS file. Important, add the final extension ".fcs".
#' @keywords FCS generation
#' @keywords final FCS
#' @keywords exporting
#' @export
#' @importFrom SummarizedExperiment colData assay
#' @examples
#' \dontrun{
#' export.metaFCS(fcs.SE = fcs_se, dr.object = dr, output.name = "final_file.fcs")
#' }

export.metaFCS <- function(fcs.SE, dr.object, output.name){
  # metadata adjusting
  mt <- sapply(colData(fcs.SE), function(x) as.numeric(as.factor(x)))
  
  # prepare dr object
  if(class(dr.object) != "list"){
    diff <- setdiff(colnames(dr.object), c(colnames(colData(fcs.SE)), rownames(fcs.SE)))
    dr <- dr.object[,colnames(dr.object) %in% diff]
  }else{
    diff <- setdiff(colnames(dr.object$dr), c(colnames(colData(fcs.SE)), rownames(fcs.SE)))
    dr <- dr.object$dr[,colnames(dr.object$dr) %in% diff]
  }
  
  # combine all together
  to_export <- cbind(t(assay(fcs.SE, i = "raw")), mt, dr)
  
  ## create FCS
  flowCore::write.FCS(premessa::as_flowFrame(as.matrix(to_export), source.frame = NULL), output.name)
}                                             


#' read.FCSset
#'
#' It reads multiple FCS files (in a computationaly-distributed way) contained in a specific folder or given through a vector.
#' @param filelist A vector with full path of FCS files to be read, commonly generated through \code{\link[base:list.files]{base::list.files()}}. If \code{NULL}, this file list will be generated as indicated below.
#' @param directory If \code{filelist = NULL}, those files stored in this location will be read. Default = \code{getwd()} (current directory).
#' @param pattern Pattern for reading files within \code{directory}. Default = \code{"fcs"}.
#' @param events Numeric vector indicating how many events to read in each FCS file (same behaviour as \code{\link[flowCore:read.FCS]{flowCore::read.FCS()}}). Default = \code{NULL}, i.e., all events will be read.
#' @param dataset Specify data segment to be read in the FCS file (same behaviour as \code{\link[flowCore:read.FCS]{flowCore::read.FCS()}}). Default = \code{1}.
#' @param num.threads Number of threads for reading the FCS files. For Windows users, you should to install \code{parallelsugar} from (\href{https:https://github.com/nathanvan/parallelsugar}{Nathanvan's GitHub}). Defult = \code{NULL}, i.e., not parallelization.
#' @keywords FCS reading
#' @keywords FCS parallel
#' @export
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
  if(events == "all") events <- NULL
  if(is.null(filelist)) print(filelist <- list.files(path = directory, pattern = pattern, full.names = T))
  
  if(Sys.info()[1] == "Windows") require(parallelsugar) #devtools::github_install('nathanvan/parallelsugar')
  fcs <- ncdfFlow::as.flowSet(ncdfFlow::read.ncdfFlowSet(filelist, emptyValue = FALSE, transformation = FALSE, truncate_max_range = FALSE, 
                                        which.lines = events, dataset = dataset, mc.cores = num.threads))
  
  return(fcs)
}


#' scMetadata.fromFCS
#'
#' It creates, from an initial metadata provided by the user, a global metadata for each event in a (\href{https://www.rdocumentation.org/packages/flowCore/versions/1.38.2/topics/flowSet-class}{flowset}) object.
#' @param flowset A flowset object generated following the (\href{https://bioconductor.org/packages/devel/bioc/vignettes/flowCore/inst/doc/HowTo-flowCore.pdf}{\code{flowCore} tutorial}).
#' @param metadata Metadata information for each read FCS file.
#' @param add.exprs Logical indicating if expression data must be added to each cell's metadata. Default = \code{TRUE}.
#' @keywords single-cell metadata
#' @keywords flowset metadata
#' @export
#' @examples
#' \dontrun{
#'  sc_metadata <- scMetadata.fromFCS(fcs, metadata, add.exprs = T)
#' }

scMetadata.fromFCS <- function(flowset, metadata, add.exprs = T){
  expr <- flowCore::fsApply(flowset, flowCore::exprs)
  mtd <- data.frame()
  pb <- progress::progress_bar$new(total = length(flowset), format = "Generating single cell metadata [:bar] :percent")
  
  file_count <- 1
  for(i in flowCore::sampleNames(flowset)){
    pb$tick()
    
    aux_md <- subset(metadata, metadata$filename == i)
    mtd_aux <- aux_md[rep(1,nrow(flowset[[i]])),]
    mtd_aux$cell_id <- paste0(file_count, ".", seq_len(nrow(mtd_aux)))
    rownames(mtd_aux) <- mtd_aux$cell_id
    
    mtd <- rbind(mtd, mtd_aux)
    file_count <- file_count + 1
    Sys.sleep(1/10)
  }
  if(add.exprs) mtd_exprs <- cbind(mtd, expr) else mtd_exprs <- mtd
  return(mtd_exprs)
}


#' exprs.saturate
#' @export

exprs.saturate <- function(expression_df){
  rng <- matrixStats::colQuantiles(expression_df, probs = c(0.01, 0.99))
  expr01 <- t((t(expression_df) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0
  expr01[expr01 > 1] <- 1
  
  return(expr01)
}


#' div.colors
#' @export
#' @importFrom RColorBrewer brewer.pal brewer.pal.info

div.colors <- function(n, set.seed = 333){
  require(RColorBrewer)
  
  if(n < 74){ #all posible colors in qual type palettes
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col = unique(col_vector)
    set.seed(set.seed); col <- sample(col, n, replace = T)
  }else{
    qual_col_pals = sample(brewer.pal.info)
    col_vector = sample(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
    col = sample(col_vector, 333)
    set.seed(set.seed); col <- sample(col, n, replace = T)
  }
  return(col)
}



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


#' fcs.SE
#'
#' It reads and creates a FCS.SE object (based on the structure of a (\href{https://www.bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html}{\code{SummarizedExperiment}}) from FCS files in a specific folder or indicated in a vector.
#' @param filelist A vector with full path of FCS files to be read, commonly generated through \code{\link[base:list.files]{base::list.files()}}. If \code{NULL}, this file list will be generated as indicated below.
#' @param directory If \code{filelist = NULL}, those files stored in this location will be read. Default = \code{getwd()} (current directory).
#' @param pattern Pattern for reading files within \code{directory}. Default = \code{"fcs"}.
#' @param events Numeric vector indicating how many events to read in each FCS file (same behaviour as \code{\link[flowCore:read.FCS]{flowCore::read.FCS()}}), those are the number of events for new and unified generated FCS files. Default = \code{NULL}, i.e., all events will be read.
#' @param dataset Specify data segment to be read in the FCS file (same behaviour as \code{\link[flowCore:read.FCS]{flowCore::read.FCS()}}). Default = \code{1}.
#' @param num.threads Number of threads for reading the FCS files. For Windows users, you should to install \code{parallelsugar} from (\href{https:https://github.com/nathanvan/parallelsugar}{Nathanvan's GitHub}). Defult = \code{NULL}, i.e., not parallelization.
#' @param metadata Metadata information for each FCS files included in the analysis. It will be used to assign to each event an specific metadata.
#' @param transformation Because FCS data is normally exported in a logarithmic format, it is necessary its transformation for downstream analysis. Default (and the only one) = \code{arcsinh}.
#' @param transf.cofactor Cofactor numeric value for that \code{arcsinh} transformation. Values can ranged from 15 to 50 (for mass cytometry) and from 150 to 10,000 (for flow cytometry). Default = \code{500}.
#' @keywords FCS.SE object
#' @keywords read FCS
#' @keywords single-cell metadata
#' @return The FCS.SE object is comprised with multiple (and downstream customizable) elements \enumerate{
#'           \item \code{colData(fcs.SE)}: A \code{data.frame} with metadata information for each event in the read FCS files.
#'           \item \code{assayNames(fcs.SE)}: Two expression matrices with markers as rows and events (cells) as columns. The first matrix (\code{raw}) contains the raw expression for each marker, and the second one (\code{transformed}), the \code{arcsinh}-transformed expression for subsequent analysis.
#'           \item \code{metadata(fcs.SE)}: Additional information about the experiment and following downstream steps: those read FCS files (\code{input_files}), \code{reduced_metada} with metadata for each FCS file...
#' @export
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
#' fcs_se <- fcs.SE(directory = "../data/", which.lines = 1000, metadata = metadata_user, 
#'     transf.cofactor = 500)
#' }

fcs.SE <- function(filelist = NULL, directory = getwd(), pattern = ".fcs$", events = "all", dataset = 1, num.threads = NULL, metadata, transformation = "arcsinh", transf.cofactor = 500){
  if(events == "all") events <- NULL

  # read FCS files
  fcs <- suppressMessages(read.FCSset(filelist, directory, pattern, events, dataset, num.threads))
  raw_data <- flowCore::fsApply(fcs, flowCore::exprs)
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
  
  # build SummarizedExperiment object
  fcs_se <- SummarizedExperiment::SummarizedExperiment(assays = list(raw = t(raw_data), transformed = t(transf_data)), colData = col_data)
  # experiment info
  fcs_se@metadata$input_fcs <- as.vector(unique(fcs_se$filename))
  fcs_se@metadata$reduced_metadata <- metadata

  return(fcs_se)  
}


#'as.flowSet.SE
#'
#' It tranforms a FCS.SE object into a (\href{https://www.rdocumentation.org/packages/flowCore/versions/1.38.2/topics/flowSet-class}{flowset}) object.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which obtain a \code{flowset} object.
#' @keywords FCS.SE to flowset
#' @keywords flowset generation
#' @export
#' @importFrom SummarizedExperiment colData assay
#' @examples
#' \dontrun{
#'  sc_metadata <- scMetadata.fromFCS(fcs, metadata, add.exprs = T)
#' }

as.flowSet.SE <- function(fcs.SE, assay.i){
  df_flowset <- list()
  for(i in unique(fcs.SE$filename)){
    aux_fcsSE <- fcs.SE[,fcs.SE$filename == i]
    df_flowset[[i]] <- premessa::as_flowFrame(t(assay(aux_fcsSE, i = "raw")))
  }
  names(df_flowset) <- fcs.SE@metadata$input_fcs
  df_flowset <- methods::as(df_flowset, "flowSet")
  return(df_flowset)
}


#'gauss.norm
#'
#' It performs a Gaussian normalization based on \code{\link[flowStats:gaussNorm]{flowStats::gaussNorm()}} function.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param marker.to.norm Marker name(s) to normalize (based on visual inspection with \code{\link[FlowCT:multidensity]{FlowCT::multidensity()}}).
#' @param norm.matrix.name New normalized matrix name (it will stored within the \code{FCS.SE} object). Default = \code{"normalized"}.
#' @keywords Gaussian normalization
#' @keywords marker alignment
#' @keywords marker normalization
#' @export
#' @examples
#' \dontrun{
#'  fcs_se <- gauss.norm(fcs.SE = fcs_se, marker.to.norm = c("CD62L", "CCR4", "SSC_A"))
#' }

gauss.norm <- function(fcs.SE, marker.to.norm, norm.matrix.name = "normalized"){
  for(i in marker.to.norm) datr <- flowStats::gaussNorm(as.flowSet.SE(fcs.SE, assay.i = "transformed"), i)$flowset
  SummarizedExperiment::assay(fcs.SE, i = norm.matrix.name) <- t(flowCore::fsApply(datr, exprs))
  return(fcs.SE)
}


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
#' fcs_se_qc <- qc.and.removeDoublets(fcs.SE = fcs_se, 
#'     physical.markers = c("FSC_A", "FSC_H", "SSC_A", "SSC_H"))
#'
#' # option 2: output an index with low quality events
#' idx_qc <- qc.and.removeDoublets(fcs.SE = fcs_se, 
#'     physical.markers = c("FSC_A", "FSC_H", "SSC_A", "SSC_H"), return.idx = T)
#' 
#' # option 3: working with all FCS files within a folder
#' idx_qc <- qc.and.removeDoublets(directory = "../data/", 
#'     physical.markers = c("FSC_A", "FSC_H", "SSC_A", "SSC_H"), output.folder = "HQ_files")
#' }

qc.and.removeDoublets <- function(fcs.SE = NULL, filelist = NULL, directory = getwd(), pattern = "fcs", physical.markers, output.folder = NULL, output.suffix = "qc", return.idx = F){
  if(!is.null(fcs.SE)){
    fcs <- as.flowSet.SE(fcs.SE, assay.i = "raw")
    filenames <- fcs@phenoData@data$name
    
    idx <- losses <- c()
    for(file in filenames){
      fcs[[file]]@description$FILENAME <- file #requirements for flow_auto_qc
      
      invisible(utils::capture.output(idx1 <- flowAI::flow_auto_qc(fcs[[file]], ChExcludeFS = physical.markers, 
        fcs_QC = FALSE, output = 3, html_report = F, folder_results = F, mini_report = F)[[1]]))
      
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


# 'clusters.rename
#'
#' It renames the numerical clusters detected through \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}} to specific cell populations identified by flow cytometry user.
#' @param x Vector with numeric values corresponding to clusters detected with \code{\link[FlowCT:fsom.clustering]{FlowCT::fsom.clustering()}}.
#' @param cluster Numeric vector or column name with numeric values for replacing.
#' @param name Cell population names to replace numeric values from \code{cluster} column.
#' @keywords population assignment
#' @keywords clusters renaming
#' @export
#' @examples
#' \dontrun{
#' # option 1  
#' head(replace_data)
#'  #   original_cluster new_cluster
#'  # 1 1                debris     
#'  # 2 2                lymphocytes
#'  # 3 3                monocytes
#'  # 4 4                lymphocytes
#'  # 5 5                eosinophils
#'  # 6 6                eosinophils
#' fcs_se$SOM_named <- clusters.rename(fcs_se$SOM, cluster = replacedata$original_cluster, 
#'     name = replacedata$new_cluster)
#' 
#' # option 2
#' fcs_se$SOM_named <- clusters.rename(fcs_se$SOM, cluster = 1:6, 
#'     name = c("debris", "lymphocytes", "monocytes", "lymphocytes", "eosinophils", "eosinophils"))
#' }

clusters.rename <- function(x, cluster, name){
  if(class(cluster) != "character") pattern <- as.character(cluster) else pattern <- cluster
  if(class(name) != "character") replacement <- as.character(name) else replacement <- name
  x2 <- as.character(x)
  
  result = x2
  for (i in 1:length(pattern)) {
    result[grep(pattern[i], x2)] = replacement[i]
  }
  if(class(x) == "factor") return(factor(result)) else return(result)
}


# 'combine.subclusterings
#'
#' It combines initial \code{FCS.SE} object (without subclustering) with other \code{FCS.SE} objects with subclustering analysis coming from downstream steps and generates a new \code{FCS.SE} object. This final \code{FCS.SE} object has an additional column combining all information from initial and subclustering analysis.
#' @param initial.fcs.SE A \code{FCS.SE} object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}. The initial one, without extracting any cell population.
#' @param subclustering.fcs.SE A list with all \code{FCS.SE} object(s) generated in the subclustering analysis (they have to come from the original \code{initial.fcs.SE}).
#' @param clusters.named Column name from the \code{initial.fcs.SE} object which contains renamed clusters (through \code{\link[FlowCT:clusters.rename]{FlowCT::clusters.rename()}}) and has been used to extract cell populations for subclustering steps.
#' @keywords final FCS.SE object
#' @keywords combine subclustering
#' @export
#' @importFrom SummarizedExperiment colData
#' @examples
#' \dontrun{
#'  fcs_se_final <- combine.subclusterings(initial.fcs.SE = fcs_se, clusters.named = "SOM_named", 
#'     subclustering.fcs.SE = list(fcs_se_lymphos, fcs_se_monos))
#' }

combine.subclusterings <- function(initial.fcs.SE, subclustering.fcs.SE, clusters.named = "SOM_named"){
  mdg <- colData(initial.fcs.SE)
  
  subclusterings <- c()
  rm.cells <- list()
  for(i in subclustering.fcs.SE){
    subclusterings <- c(subclusterings, i@metadata$subclustering)
    
    md_sub <- colData(i)
    
    # add differential cols from subclusterings
    diff1 <- setdiff(colnames(md_sub), colnames(mdg))
    mdg[,diff1] <- NA
    diff2 <- setdiff(colnames(mdg), colnames(md_sub))
    md_sub[,diff2] <- NA
    
    # extract those subclustered samples from original fcs.SE object
    mdg <- rbind(mdg[mdg[,clusters.named] != i@metadata$subclustering,], md_sub)
    
    # delete removed populations from the (i)th subclustering within the first fcs.SE object
    if(!is.null(i@metadata$removed_populations)){
      for(j in names(i@metadata$removed_populations))
        rm.cells[[j]] <- i@metadata$removed_populations[[j]]
        mdg <- mdg[setdiff(mdg$cell_id, i@metadata$removed_populations[[j]]),]
        
        # substract deleted cells from original assays
        initial.fcs.SE <- initial.fcs.SE[, initial.fcs.SE$cell_id %in% setdiff(mdg$cell_id, i@metadata$removed_populations[[j]])]
    }
  }
  
  # create final named_cluster col ---> beta, this will bring problems with multiple subclusterings!
  diff <- setdiff(colnames(mdg), colnames(colData(initial.fcs.SE)))
  for(i in diff){
    if(suppressWarnings(sum(is.na(as.numeric(as.character(mdg[,i])))) == nrow(mdg))){ #detect names_clusters cols in subclustering
      mdg$tmp <- ifelse(is.na(mdg[,i]), as.character(mdg[,clusters.named]), 
                        as.character(mdg[,i]))
    }
    mdg[,i] <- ifelse(is.na(mdg[,i]), 0, mdg[,i]) #replace NAs by 0 to avoid FCS wrong building
  }
  colnames(mdg)[ncol(mdg)] <- paste0(clusters.named, "_final")
  
  SummarizedExperiment::colData(initial.fcs.SE) <- mdg
  initial.fcs.SE@metadata$subclusterings$populations <- paste(subclusterings, collapse = " + ")
  initial.fcs.SE@metadata$subclusterings$removed_populations <- rm.cells
  
  return(initial.fcs.SE)
}


# 'remove.pop
#'
#' It removes one or multiple cell populations from a \code{FCS.SE} object.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param population Name(s) of cell population(s) to be removed.
#' @param clusters.named Column name from the \code{colData(fcs.SE)} object which contains renamed clusters (through \code{\link[FlowCT:clusters.rename]{FlowCT::clusters.rename()}}).
#' @keywords remove cell population
#' @keywords debris
#' @export
#' @importFrom SummarizedExperiment colData
#' @examples
#' \dontrun{
#' fcs_se_rm <- remove.pop(fcs_se, clusters.named = "SOM_named", 
#'     population = c("debris", "unclassified"))
#' }

remove.pop <- function(fcs.SE, population, clusters.named){
  for(i in population){
    fcs.SE@metadata$removed_populations[[i]] <- colData(fcs.SE)[fcs.SE[[clusters.named]] == i,"cell_id"]
    fcs.SE <- fcs.SE[,fcs.SE[[clusters.named]] != i]
    fcs.SE[[clusters.named]] <- droplevels(fcs.SE[[clusters.named]])
  }
  return(fcs.SE)
}


# 'median.values
#'
#' It calculates median values according a specfied variable, normally "filename".
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate medians. Default = \code{"raw"}.
#' @param var Variable for grouping and calculating medians. Default = \code{"filename"} (i.e., names of FCS files).
#' @keywords median values
#' @keywords MFI
#' @keywords median fluorescence intensity
#' @export median.values
#' @examples
#' \dontrun{
#' mfis_FCS <- median.values(fcs.SE = fcs_se)
#' med_SOM_clust <- median.values(fcs.SE = fcs_se, assay.i = "normalized", var = "SOM_named")
#' }

median.values <- function(fcs.SE, assay.i = "raw", var = "filename"){
  med <- list()
  for(i in unique(fcs.SE[[var]])){
    aux_se <- fcs.SE[,fcs.SE[[var]] == i]
    med[[i]] <- apply(SummarizedExperiment::assay(aux_se, i = assay.i), 1, stats::median)
  }
  if(is.null(names(med))) names(med) <- unique(fcs.SE[[var]])
  return(t(as.data.frame(med)))
}

# 'marker.names
#'
#' It shows marker names of a \code{FCS.SE} object and renames them accordin a new vector provided by the user.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param new.names Vector with new channel/marker names (it must has the same length that \code{FCS.SE}'s markers). Default = \code{NULL} (i.e., markers will not be renamed, only displayec).
#' @keywords median values
#' @keywords MFI
#' @keywords median fluorescence intensity
#' @export median.values
#' @examples
#' \dontrun{
#' marker.names(fcs_se)
#' new_names <- c("FSC_A", "FSC_H", "SSC_A", "SSC_H", "CD62L", "CXCR3", "CD8", "CCR4", "CCR6", "CD4", "CD45", "CD27")
#' fcs_se <- marker.names(fcs_se, new.names = new_names)
#' }

marker.names <- function (fcs.SE, new.names = NULL){
  raw_names <- rownames(fcs.SE)
  if(is.null(new.names)){
    if(is.null(fcs.SE@metadata$raw_channel_names)){
      print(knitr::kable(data.frame(raw_names, "not defined yet!"), col.names = c("raw name", "new name")))
    }else{
      print(knitr::kable(data.frame(fcs.SE@metadata$raw_channel_names, raw_names), 
                         col.names = c("raw name", "new name")))
    }
  }else{
    ## add new markers names
    if(length(new.names) != length(raw_names)) stop("New names must to have the same length to the original ones!")
    print(knitr::kable(data.frame(raw_names, new.names), col.names = c("raw name", "new name")))
    
    rownames(fcs.SE) <- new.names
    fcs.SE@metadata$raw_channel_names <- raw_names
    return(fcs.SE)
  }
}

# ### Accessor methods ###############################################################################
# #' @exportMethod metadata
# setMethod("metadata", "FlowCT", function(x)
#     slot(x, "metadata"))

# #' @exportMethod colData
# setMethod("colData", "FlowCT", function(x, ...) {
#     slot(x, "colData")
# })
