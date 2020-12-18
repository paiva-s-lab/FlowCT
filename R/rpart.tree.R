#' Recursive partitioning tree
#'
#' This function performs a recursive partitioning tree using the percentage (or raw counts) of identified cell populations. Depending if \code{time.var} is present or absent, final tree were a cassification or regression tree, respectively.
#' @param fcs.SCE A \code{fcs.SCE} object generated through \code{\link[FlowCT:fcs.SCE]{FlowCT::fcs.SCE()}}.
#' @param cell.clusters Name of column containing clusters identified through \code{\link[FlowCT:clustering.flow]{FlowCT::clustering.flow()}}.
#' @param variables Vector with variables for calculating the cutoff. If nothing is detailed (\code{NULL}, default), all immune populations from \code{cell.clusters} will be considered.
#' @param value String specifying if final resuls should be proportions ("percentage", default) or raw counts ("counts").
#' @param time.var Survival time variable.
#' @param event.var Variable with event censoring.
#' @param xerror Numeric value with cross-validation error for pruning the tree. If missing, the minimal error will be automatically selected.
#' @param return.data Original, pruned and used metadata should be returned?. Default = \code{FALSE}.
#' @keywords recursive partioning
#' @keywords CART
#' @keywords classification regression tree
#' @export rpart.tree
#' @import dplyr
#' @importFrom SingleCellExperiment colData
#' @importFrom stats as.formula
#' @importFrom graphics par
#' @examples
#' \dontrun{
#' tr <- rpart.tree(fcs.SCE = fcs, cell.clusters = "clusters_named", 
#'    time.var = "PFS", event.var = "PFS_c", return.data = t)
#' }


rpart.tree <- function(fcs.SCE, cell.clusters, variables, value = "percentage", time.var, event.var, xerror, return.data = F){
  if (!requireNamespace(c("rpart", "rpart.plot"), quietly = TRUE)) stop("Packages \"rpart\" and \"rpart.plot\" needed for this function to work. Please install them.", call. = FALSE)
  
  cell_props <- as.data.frame.matrix(t(barplot.cell.pops(fcs.SCE, cell.clusters = cell.clusters, plot = F, count.by = "filename", return.mode = value)))
  md <- colData(fcs.SCE) %>% as.data.frame() %>% distinct(.data$filename, .keep_all = T)
  surv_data <- data.frame(cell_props[match(md$filename, rownames(cell_props)),], md)

  if(missing(variables)) variables <- fcs.SCE[[cell.clusters]]
  if(class(surv_data[,event.var]) != "numeric") surv_data[,event.var] <- as.numeric(as.factor(surv_data[,event.var]))

  set.seed(33)
  if(missing(time.var)){
    f <- as.formula(paste0(event.var, " ~ ", paste(unique(variables), collapse = "+")))
    tr <- rpart::rpart(f, data = surv_data, method = "class", control = rpart::rpart.control(cp = 0))
  }else{
    if(class(surv_data[,time.var]) != "numeric") surv_data[,time.var] <- as.numeric(surv_data[,time.var])
    
    cat(">>> Attention,", as.numeric(table(surv_data[,time.var] == 0)[2]), "value(s) were discarded because", time.var, "was zero.\n")
    surv_data <- surv_data[surv_data[,time.var] > 0,]
    
    f <- as.formula(paste0("Surv(", time.var, ", ", event.var, ") ~ ", paste(unique(variables), collapse = "+")))
    tr <- rpart::rpart(f, data = surv_data, method = "exp", control = rpart::rpart.control(cp = 0))
  }

  if(missing(xerror)){ #minimal crossvalidated error
    xerror <- tr$cptable[which.min(tr$cptable[,"xerror"]),"CP"]
    tr_p <- rpart::prune(tr, cp = xerror)
  }else{
    tr_p <- rpart::prune(tr, cp = xerror)
  }

  par(mfrow = c(1,2)) #combine plots
  rpart.plot::rpart.plot(tr, under = T, type = 5, extra = 102, nn = F, cex = 0.7, main = "original")
  rpart.plot::rpart.plot(tr_p, under = T, type = 5, nn = F, extra = 102, cex = 0.7, main = paste0("pruned (", round(xerror, 3), ")")) 

  if(return.data) return(list(original_tree = tr, pruned_tree = tr_p, metadata = md))
}
