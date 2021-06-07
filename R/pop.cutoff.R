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
#' @export extract.cutoffs
#' @import dplyr
#' @importFrom stats quantile median
#' @examples
#' \dontrun{
#' ct <- pop.cutoff(fcs.SCE = fcs, cell.clusters = "SOM_named", time.var = "PFS",
#'     event.var = "PFS_c", cutoff.type = "quantiles")
#'     
#' extract.cutoffs(ct)     
#' }

prog.pop.selection <- function(fcs.SCE, assay.i = "normalized", cell.clusters, variables, cutoff.type = "maxstat", time.var, event.var, condition.col, cell.value = "percentage", method, method.params, plot = T, return.ML.object = F, train.index){
  require(survival)
  
  if(class(fcs.SCE[[time.var]]) != "numeric" || class(fcs.SCE[[event.var]]) != "numeric") 
    stop(call. = F, 
         "Please, your time and event variables must be in numeric format... IMPORTANT: positive or negative event coded as 1 and 0, respectively.")
  
  if(missing(variables)){
    variables <- levels(unique(fcs.SCE[[cell.clusters]]))
    suppressWarnings(if(!is.na(sum(as.numeric(variables)))) variables <- paste0("X", variables) else variables <- variables) #in case variable is not named yet (only cluster number)
  }
  
  if(missing(condition.col)) condition.col <- NULL
  
  ## train/test datasets --->>> TODO: extremelly bigger loop, simplify?
  if(!missing(train.index)){
    fcs <- list(train = fcs.SCE[,fcs.SCE$filename %in% train.index], 
                test = fcs.SCE[,!(fcs.SCE$filename %in% train.index)])
    
    ## cutoff calculation
    if(cutoff.type != "none" && cutoff.type %in% c("maxstat", "median", "quantiles", "quantile", "terciles", "tercile")){
      dataset_surv <- lapply(fcs, function(x) 
        pop.cutoff(fcs.SCE = x, cell.clusters = cell.clusters, time.var = time.var, event.var = event.var, value = cell.value, cutoff.type = cutoff.type, assay.i = assay.i))
      cts <- lapply(dataset_surv, function(x) extract.cutoffs(x))
      
      dataset_surv <- lapply(dataset_surv, function(x) x[!(colnames(x) %in% variables)])
      variables <- paste0(variables, ".c")
      
      
    }else if(cutoff.type == "none"){
      prop_table_surv <- lapply(fcs, function(x) 
        barplot.cell.pops(fcs.SCE = x, cell.clusters = cell.clusters, count.by = "filename", return.mode = cell.value, plot = F, assay.i = "normalized"))
      
      dataset_surv <- lapply(prop_table_surv, function(x) 
        merge(dplyr::distinct(as.data.frame(colData(fcs.SCE)), .data$filename, .keep_all = T), 
              data.frame(filename = colnames(x), as.data.frame.matrix(t(x))), 
              by = "filename"))
      
    }else{
      stop(call. = F, "Please, select a valid method (see ?pop.cutoff): 'mean', 'quantiles' or 'maxstat' (default), or 'none' for not cutoff.type calculation.")
    }
    
    ## recast to numeric
    dataset_surv <- lapply(dataset_surv, function(x){
      num <- as.data.frame(apply(x[,c(time.var, event.var, grep(paste0("^", variables, "$", collapse = "|"), colnames(x), value = T))], 2, 
                                 function(y) if(class(y) != "numeric") as.numeric(as.factor(y)) else y))
      rownames(num) <- rownames(x)
      num$condition <- unlist(x[condition.col]) #condition.col must be character/factor for biosigner
      return(num)
    })
    
  }else{ #>> if no train/test
    ## cutoff calculation
    if(cutoff.type != "none" && cutoff.type %in% c("maxstat", "median", "quantiles", "quantile", "terciles", "tercile", "roc")){
      dataset_surv <- pop.cutoff(fcs.SCE = fcs.SCE, cell.clusters = cell.clusters, time.var = time.var, event.var = event.var, value = cell.value, cutoff.type = cutoff.type, assay.i = assay.i)
      cts <- extract.cutoffs(dataset_surv)
      
      dataset_surv <- dataset_surv[,!(colnames(dataset_surv) %in% variables)]
      variables <- paste0(variables, ".c")
      
      
    }else if(cutoff.type == "none"){
      prop_table_surv <- barplot.cell.pops(fcs.SCE = fcs.SCE, cell.clusters = cell.clusters, count.by = "filename", return.mode = cell.value, plot = F, assay.i = "normalized")
      
      dataset_surv <- merge(dplyr::distinct(as.data.frame(colData(fcs.SCE)), .data$filename, .keep_all = T), 
                            data.frame(filename = colnames(prop_table_surv), as.data.frame.matrix(t(prop_table_surv))), 
                            by = "filename")
      
    }else{
      stop(call. = F, "Please, select a valid method (see ?pop.cutoff): 'mean', 'quantiles' or 'maxstat' (default), or 'none' for not cutoff.type calculation.")
    }
    
    ## recast to numeric
    num <- as.data.frame(apply(dataset_surv[,c(time.var, event.var, grep(paste0("^", variables, "$", collapse = "|"), colnames(dataset_surv), value = T))], 2,
                               function(y) if(class(y) != "numeric") as.numeric(as.factor(y)) else y))
    rownames(num) <- dataset_surv$filename
    num$condition <- unlist(dataset_surv[condition.col]) #condition.col must be character/factor for biosigner
    dataset_surv <- list(train = num) #coerce to list for ML downstream (built for train/test)
  }
   
  ## ML functions  
  if(method == "biosign"){
    if (!requireNamespace("biosigner", quietly = TRUE)) stop("Package \"biosigner\" is needed for this function. Please install it.", call. = FALSE)
    
    if(!missing(method.params)) 
      res <- do.call(biosigner::biosign, c(list(x = dataset_surv$train[grep(paste0("^", variables, "$", collapse = "|"), colnames(dataset_surv$train))], 
                                                y = as.vector(unlist(dataset_surv$train[condition.col]))), method.params)) 
    else
      res <- biosigner::biosign(x = dataset_surv$train[grep(paste0("^", variables, "$", collapse = "|"), colnames(dataset_surv$train))], 
                                y = as.vector(unlist(dataset_surv$train[condition.col])))
    
    features_selection <- res@signatureLs[-4]
    
  }else if(method == "random_forest"){
    if (!requireNamespace("randomForestSRC", quietly = TRUE)) stop("Package \"randomForestSRC\" is needed for this function. Please install it.", call. = FALSE)
    
    f <- as.formula(paste0("Surv(", time.var, ",", event.var, ") ~ ", 
                           paste(grep(paste0("^", variables, "$", collapse = "|"), colnames(dataset_surv$train), value = T), collapse = " + ")))
    
    if(!missing(method.params)) 
      res <- do.call(randomForestSRC::rfsrc, c(list(formula = f, data = dataset_surv$train, seed = 333), method.params)) 
    else
      res <- randomForestSRC::rfsrc(f, data = dataset_surv$train, seed = 333)
    
    features_selection <- merge(as.data.frame(randomForestSRC::vimp(res)$importance), 
                         randomForestSRC::var.select(res, verbose = F)$varselect[,1, drop = F], 
                         by = "row.names")
    colnames(features_selection) <- c("variable", "VIMP", "depth")
    
    if(plot){
      print(ggplot(features_selection, aes(VIMP, reorder(variable, VIMP), fill = depth)) + 
        geom_bar(stat = "identity") + 
        ylab("Cell population") + xlab("(-) event related <--- VIMP ---> (+) event related") + 
        theme_bw())
    }
    
  }else if(tolower(method) == "survboost"){ ###IMPORTANT: for using SurvBoost, you MUST to library it BEFORE FlowCT
    if (!requireNamespace("SurvBoost", quietly = TRUE)) stop("Package \"SurvBoost\" is needed for this function. Please install it ---> devtools::install_github(\"EmilyLMorris/survBoost\"", call. = FALSE)
    
    f <- as.formula(paste0("Surv(", time.var, ",", event.var, ") ~ ", paste(grep(paste0("^", variables, "$", collapse = "|"), colnames(dataset_surv$train), value = T), collapse = " + ")))
    
    if(!missing(method.params)) 
      res <- rlang::invoke(SurvBoost::boosting_core, c(list(formula = f, data = as.data.frame(dataset_surv$train)), method.params)) 
    else
      res <- SurvBoost::boosting_core(f, data = as.data.frame(dataset_surv$train), rate = 0.1)
    
    features_selection <- res$coefficients[res$coefficients != 0]
    
    if(plot){
      selection_df <- rbind(rep(0, ncol(res$selection_df)), res$selection_df)
      colnames(selection_df) <- names(res$coefficients)
      plot_data <- reshape::melt(data.frame(x_axis = c(0:res$mstop), selection_df), id.vars = c("x_axis"))
      
      print(ggplot(plot_data, aes(x = x_axis, y = value, group = variable)) + geom_line() + theme_bw() + 
              theme(text = element_text(size = 16)) + 
              ylab("Coefficient Estimate") + xlab("Number of iterations") + 
              ggrepel::geom_label_repel(data = plot_data[plot_data$x_axis == res$mstop,], aes(label = variable, x = x_axis, y = value), max.overlaps = length(res$coefficients)))
    }
    
  }else{
    stop("Please, indicate a valid method: 'biosign', 'random_forest' or 'survboost'.", call. = F)
  }
  
  ## returning
  if(return.ML.object){
    if(missing(method.params)) 
      return(list(ML.object = res, 
                  survival.data = list(data = dataset_surv, cutoffs = cts))) else
      return(list(ML.object = res, 
                  survival.data = list(data = dataset_surv, cutoffs = cts), 
                  method.params = data.frame(method.params = unlist(method.params))))
  }else{
    if(missing(method.params)) 
      return(features_selection) else
      return(list(features_selection = features_selection, 
                  method.params = data.frame(method.params = unlist(method.params))))
  }
}

                               
### class and methods
setClass("cutoff.object", slots = list(data.c = "data.frame", cutoffs = "vector"))

setMethod("show", "cutoff.object", function(x) print(x@data.c))
setMethod("colnames", "cutoff.object", function(x) colnames(x@data.c))
setMethod("head", "cutoff.object", function(x) head(x@data.c))
setMethod("tail", "cutoff.object", function(x) tail(x@data.c))
setMethod("[", c("cutoff.object", "ANY"), function(x,i,..., drop = T) x@data.c = x@data.c[i])

#' @name extract.cutoffs
#' @docType methods
#' @rdname extract-methods
#'
setGeneric("extract.cutoffs", function(x) standardGeneric("extract.cutoffs"))
setMethod("extract.cutoffs", "cutoff.object", function(x) x@cutoffs)
