#' @title 
#' Determine optimal cutoffs for immune populations
#'
#' @description 
#' Once immune cell populations are identified, this function calculates a cutoff from the percentage (or raw counts).
#' 
#' @details
#' There four available methods for cytoff calculation:\itemize{
#' \item \code{median}
#' \item \code{quantiles} is leveraged on quantile categorization.
#' \item \code{maxstat} (default), based on \href{https://rpkgs.datanovia.com/survminer/reference/surv_cutpoint.html}{maximally selected ranks statistics}. It takes into account the survival time and censoring event for cutoff calculation.
#' \item \code{roc}, the classical ROC-based calculation according Youden's index. It considerates the censoring event for categorizing.
#' }
#' 
#' Note: for accessing to cutoff values used for categorization...
#' \code{`extract.cutoffs(your_cutoff_object)`}
#' 
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param cell.clusters Name of column containing clusters identified through \code{\link[FlowCT:clustering.flow]{FlowCT::clustering.flow()}}.
#' @param value String specifying if final resuls should be proportions ("percentage", default) or raw counts ("counts").
#' @param time.var Survival time variable.
#' @param event.var Variable with event censoring.
#' @param variables Vector with variables for calculating the cutoff. If nothing is detailed (\code{NULL}, default), all immune populations from \code{cell.clusters} will be considered.
#' @param cutoff.type Method for calculating survival cutoffs. Available methods are "maxstat" (default){\code{maxstat}}), "ROC", "quantiles" (i.e., terciles) and "median". 
#' @keywords survival cutoffs
#' @export pop.cutoff
#' @import dplyr
#' @importFrom stats quantile median
#' @examples
#' \dontrun{
#' ct <- pop.cutoff(fcs.SCE = fcs, cell.clusters = "SOM_named", time.var = "PFS",
#'     event.var = "PFS_c", cutoff.type = "quantiles")
#'     
#' extract.cutoffs(ct)     
#' }

pop.cutoff <- function(fcs.SCE, assay.i = "normalized", cell.clusters, value = "percentage", time.var, event.var,
                       cutoff.type = "maxstat", variables){
  ## prepare data
  prop_table_surv <- barplot.cell.pops(fcs.SCE, cell.clusters, count.by = "filename", return.mode = value, plot = F, assay.i = assay.i)
  dataset_surv <- merge(distinct(as.data.frame(colData(fcs.SCE)), .data$filename, .keep_all = T), 
                        as.data.frame.matrix(t(prop_table_surv)), by.x = "filename", by.y = "row.names")
  
  if(missing(variables)) groups <- as.character(unique(fcs.SCE[[cell.clusters]])) else groups <- variables
  
  ## cutoffs
  if(tolower(cutoff.type) == "maxstat"){
    if(class(dataset_surv[,time.var]) != "numeric") dataset_surv[,time.var] <- as.numeric(dataset_surv[,time.var])
    if(class(dataset_surv[,event.var]) != "numeric") dataset_surv[,event.var] <- as.numeric(as.factor(dataset_surv[,event.var]))
    
    if (!requireNamespace("survminer", quietly = TRUE)) stop("Package \"survminer\" is needed for this function. Please install it.", call. = FALSE)
    
    res_cut <- survminer::surv_cutpoint(dataset_surv, time = time.var, event = event.var, variables = groups, minprop = 0.2, progressbar = F)
    res_cat <- data.frame(dataset_surv, survminer::surv_categorize(res_cut)[,-c(1:2)]) #not to duplicate time and event cols
    colnames(res_cat)[grep("\\.1", colnames(res_cat))] <- gsub("\\.1", ".c", colnames(res_cat[grep("\\.1", colnames(res_cat))]))
    
    cts <- setNames(res_cut$cutpoint$cutpoint, rownames(res_cut$cutpoint))
    
  }else if(tolower(cutoff.type) == "roc"){
    if (!requireNamespace("cutpointr", quietly = TRUE)) stop("Package \"cutpointr\" is needed for this function. Please install it.", call. = FALSE)
    
    if(class(dataset_surv[,event.var]) != "numeric") dataset_surv[,event.var] <- as.numeric(as.factor(dataset_surv[,event.var]))
    
    ct <- cutpointr::multi_cutpointr(data = dataset_surv, x = groups, class = !!event.var, #shorturl.at/jC057
                                     na.rm = T, metric = cutpointr::youden, silent = T)
    
    res_cat <- sapply(1:length(groups), function(x) factor(ifelse(dataset_surv[,groups[x]] >= ct$optimal_cutpoint[x], "high", "low")))
    colnames(res_cat) <- paste0(groups, ".c")
    res_cat <- data.frame(dataset_surv, res_cat)
    
    cts <- setNames(ct$optimal_cutpoint, ct$predictor)
    
  }else if(cutoff.type %in% c("quantiles", "quantile", "tercile", "terciles")){
    res_cat <- data.frame(dataset_surv, sapply(dataset_surv[,groups], function(x){
      aux <- quantile(x, probs = seq(0,1,0.333))
      factor(ifelse(x > aux[3], "high", ifelse(x < aux[2], "low", "mid")))}))
    colnames(res_cat)[grep("\\.1", colnames(res_cat))] <- gsub("\\.1", ".c", colnames(res_cat[grep("\\.1", colnames(res_cat))]))
    
    cts <- sapply(dataset_surv[,groups], function(x) paste(round(quantile(x, probs = seq(0,1,0.333))[2:3], 2), collapse = " - "))
    
  }else if(tolower(cutoff.type) == "median"){
    res_cat <- data.frame(dataset_surv, sapply(dataset_surv[,groups], function(x) factor(ifelse(x >= median(x), "high", "low"))))
    colnames(res_cat)[grep("\\.1", colnames(res_cat))] <- gsub("\\.1", ".c", colnames(res_cat[grep("\\.1", colnames(res_cat))]))
    
    cts <- sapply(dataset_surv[,groups], function(x) round(median(x), 2))
    
  }else{
    stop("Please, specify one valid option: maxstat, ROC, quantiles or median.", call. = F)
  }
  
  return(new("cutoff.object", data.c = res_cat, cutoffs = cts))
}


### class and methods
setClass("cutoff.object", slots = list(data.c = "data.frame", cutoffs = "vector"))

setMethod("show", "cutoff.object", function(x) print(x@data.c))
setMethod("colnames", "cutoff.object", function(x) colnames(x@data.c))
setMethod("head", "cutoff.object", function(x) head(x@data.c))
setMethod("tail", "cutoff.object", function(x) tail(x@data.c))
setMethod("[", c("cutoff.object", "ANY"), function(x,i,..., drop = T) x@data.c = x@data.c[i])

setGeneric("extract.cutoffs", function(x) standardGeneric("extract.cutoffs"))
setMethod("extract.cutoffs", "cutoff.object", function(x) x@cutoffs)
