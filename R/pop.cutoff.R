#' Immune cutoffs
#'
#' Once immune cell populations are identified, it calculates a cutoff from the percentage (or raw counts) according survival time.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param cell.clusters Name of column containing clusters identified through \code{\link[FlowCT:clustering.flow]{FlowCT::clustering.flow()}}.
#' @param value String specifying if final resuls should be proportions ("percentage", default) or raw counts ("counts").
#' @param time.var Survival time variable.
#' @param event.var Variable with event censoring.
#' @param variables Vector with variables for calculating the cutoff. If nothing is detailed (\code{NULL}, default), all immune populations from \code{cell.clusters} will be considered.
#' @param cutoff.type Method for calculating survival cutoffs. Available methods are "maxstat" (default, based on \href{https://cran.r-project.org/web/packages/maxstat/index.html}{\code{maxstat}}), "quantiles" (i.e., terciles) and "median". 
#' @keywords survival cutoffs
#' @export pop.cutoff
#' @import dplyr
#' @importFrom stats quantile median
#' @examples
#' \dontrun{
#' ct <- pop.cutoff(fcs.SCE = fcs, cell.clusters = "SOM_named", time.var = "PFS",
#'     event.var = "PFS_c", cutoff.type = "quantiles")
#' }

pop.cutoff <- function(fcs.SCE, assay.i = "normalized", cell.clusters, value = "percentage", time.var, event.var, 
                       cutoff.type = "maxstat", variables){
  ## prepare data
  prop_table_surv <- barplot.cell.pops(fcs.SCE, cell.clusters, count.by = "filename", return.mode = value, plot = F, assay.i = assay.i)
  dataset_surv <- merge(distinct(as.data.frame(colData(fcs.SCE)), .data$filename, .keep_all = T), 
                        as.data.frame.matrix(t(prop_table_surv)), by.x = "filename", by.y = "row.names")
  
  if(class(dataset_surv[,time.var]) != "numeric") dataset_surv[,time.var] <- as.numeric(dataset_surv[,time.var])
  if(class(dataset_surv[,event.var]) != "numeric") dataset_surv[,event.var] <- as.numeric(as.factor(dataset_surv[,event.var]))
  if(missing(variables)) groups <- as.character(unique(fcs.SCE[[cell.clusters]])) else groups <- variables
  
  ## cutoffs
  if(cutoff.type == "maxstat"){
    if (!requireNamespace("survminer", quietly = TRUE)) stop("Package \"survminer\" needed for this function to work. Please install it.", call. = FALSE)

    res_cut <- survminer::surv_cutpoint(dataset_surv, time = time.var, event = event.var, variables = groups, minprop = 0.2, progressbar = F)
    res_cat <- data.frame(dataset_surv, survminer::surv_categorize(res_cut)[,-c(1:2)]) #not to duplicate time and event cols
    colnames(res_cat)[grep("\\.1", colnames(res_cat))] <- paste0(gsub("\\.1", "", colnames(res_cat[grep("\\.1", colnames(res_cat))])), "..", round(res_cut$cutpoint$cutpoint, 2))
  }else if(cutoff.type == "quantiles"){
    res_cat <- data.frame(dataset_surv, sapply(dataset_surv[,groups], function(x){
      aux <- quantile(x, probs = seq(0,1,0.333))
      factor(ifelse(x > aux[3], "high", ifelse(x < aux[2], "low", "mid")))}))

    cts <- sapply(dataset_surv[,groups], function(x) paste(round(quantile(x, probs = seq(0,1,0.333))[2:3], 2), collapse = "_"))
    colnames(res_cat)[grep("\\.1", colnames(res_cat))] <- paste0(gsub("\\.1", "", colnames(res_cat[grep("\\.1", colnames(res_cat))])), "..", cts)
    
    res_cat[,grepl("\\.[0-9]*_[0-9]*", colnames(res_cat))] <- lapply(res_cat[,grepl("\\.[0-9]*_[0-9]*", colnames(res_cat))], function(x) factor(x, levels = c("high", "mid", "low")))
  }else if(cutoff.type == "median"){
    res_cat <- data.frame(dataset_surv, sapply(dataset_surv[,groups], function(x) factor(ifelse(x >= median(x), "high", "low"))))

    cts <- sapply(dataset_surv[,groups], function(x) round(median(x), 2))
    colnames(res_cat)[grep("\\.1", colnames(res_cat))] <- paste0(gsub("\\.1", "", colnames(res_cat[grep("\\.1", colnames(res_cat))])), "..", cts)
   }else{
    stop("Please, specify one valid option: maxstat, quantiles or median.", call. = F)
  }
  
  return(res_cat)
}
