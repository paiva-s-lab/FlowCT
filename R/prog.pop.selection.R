#' @title 
#' Determine populations with prognostic value
#'
#' @description 
#' Determine, by using different machine learning approaches, those populations with prognostic value... either with raw percentage (continuous variable) or under cutoff.types (categorical variable).
#' 
#' @details 
#' Up to now, this wrapper function is comprising three different methods. Please, check each package's help for further details.\itemize{
#' \item \href{http://www.bioconductor.org/packages/release/bioc/vignettes/biosigner/inst/doc/biosigner-vignette.html}{\pkg{biosigner}}. It includes three classification algorithms: PLS-DA, RF and SVM; it only works with censoring event variable (but not with survival time).
#' \item \href{https://kogalur.github.io/randomForestSRC/theory.html}{\pkg{randomForestSRC}}. Random Forests for survival (and regression and classification).
#' \item \href{https://github.com/EmilyLMorris/survBoost}{\pkg{SurvBoost}}. A high dimensional variable selection method for stratified proportional hazards model.
#' }
#' 
#' The returning object changes according chosen arguments:\itemize{
#' \item if \code{return.ML.object = FALSE}, only variables' importance/coefficients will be showed;
#' \item if \code{return.ML.object = TRUE} and \code{train.index} is empty, the object (list-type) includes the machine learning object and a double data.frame with variables' importance/coefficients; and 
#' \item if \code{return.ML.object = TRUE} and \code{train.index} contains a vector of samples, the list-type object would store also a double data.frame with training and validation (test) datasets (with percentage or categorized data).
#' }
#' 
#' \strong{Important note}: for using SurvBoost's method, you MUST to load it BEFORE FlowCT for avoiding internal conflicts.
#' 
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param assay.i Name of matrix stored in the \code{fcs.SCE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param cell.clusters Name of column containing clusters identified through \code{\link[FlowCT:clustering.flow]{FlowCT::clustering.flow()}}.
#' @param variables Vector with variables for calculating the prognostic relevance. If nothing is detailed (default), all immune populations from \code{cell.clusters} will be considered.
#' @param cutoff.type Method for calculating survival cutoff.types. Available methods are "maxstat" (default), "ROC", "quantiles" (i.e., terciles) and "median". If "none" is especified, raw percentages (or counts) were used instead of categorical variables.
#' @param time.var Survival time variable.
#' @param event.var Variable with event censoring. Important note: positive and negative events should be coded as 1 and 0, respectively.
#' @param condition.col Variable with differential condition (only needed if \code{method = "biosign"}).
#' @param cell.value String specifying if final resuls should be proportions ("percentage", default) or raw counts ("counts").
#' @param method Machine learning approaches available for variable selection. Possible values are: "biosign", "random_forest" and "survboost".
#' @param method.params Internal options for "tunning" the selected method, see each package's help for more information and default values. 
#' @param plot Whether results should be plotted. Default = \code{TRUE}.
#' @param return.ML.object Logical indicating if the machine learning object must be returned (for later \code{predict}s). Default = \code{FALSE}.
#' @param train.index Vector (based on "filename" variable) with samples selected as training dataset (needed for later \code{predict}s).
#' @keywords survival prognostic
#' @export prog.pop.selection
#' @import survival
#' @import dplyr
#' @examples
#' \dontrun{
#' # eg1: only return more implied populations, after cutoff calculation
#' ml1 <- prog.pop.selection(fcs.SCE = fcs, cell.clusters = "SOM_named", 
#'           time.var = "PFS", event.var = "PFS_c", cutoff.type = "quantiles", 
#'           method = "survboost", method.params = list(rate = 0.4))
#'     
#' # eg2: apply predict (with training and validation datasets, 70%/30%), no cutoffs
#' train_idx <- sample(length(fcs$patient_id), length(fcs$patient_id)*0.7)
#' ml2 <- prog.pop.selection(fcs.SCE = fcs, cell.clusters = "SOM_named", 
#'           time.var = "PFS", event.var = "PFS_c", cutoff.type = "none", 
#'           method = "random_forest", train.index = train_idx, return.ML.object = T)
#' ml2_pr <- predict(ml2$ML.object, newdata = ml2$survival.data$test)
#' 
#' biosigner::predict(ml2$ML.object, #predict for biosigner
#'            newdata = ml2$survival.data$test[,-c(1:2,32)]) #delete survival and condition cols
#'            
#' SurvBoost::predict.boosting(ml2$ML.object, newdata = ml2$survival.data$test) #survboost           
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
