#' Survival plots according cell cutoffs
#'
#' Once immune cell cutoffs are calculated, this function draw the subsequent Kaplan-Meier curves for specified populations.
#' @param pop.cutoff.obj An object generated through \code{\link[FlowCT:pop.cutoff]{FlowCT::pop.cutoff()}}.
#' @param time.var Survival time variable.
#' @param event.var Variable with event censoring.
#' @param palette Colors vector. Default = \code{"jco"}.
#' @param variables Vector with variables for drawing the survival curve. If nothing is detailed (\code{NULL}, default), all cutoffs will be plotted.
#' @param curve.type Plotting curve methodology, i.e. "survival" (default) or "cumulative" (for progressions). 
#' @keywords survival cutoffs Kaplan-Meier curves
#' @export surv.plot
#' @importFrom stats as.formula
#' @examples
#' \dontrun{
#' surv.plot(pop_cuts, time.var = "PFS", event.var = "PFS_c", curve.type = "survival")
#' }


surv.plot <- function(pop.cutoff.obj, time.var, event.var, palette = "jco", variables, curve.type = "survival"){
  if (!requireNamespace(c("survminer", "ggpubr", "cowplot"), quietly = TRUE)) stop("Packages \"survminer\", \"cowplot\" and \"ggpubr\" needed for this function to work. Please install them.", call. = FALSE)

  if(missing(variables)) variables <- grep(".ct", colnames(pop.cutoff.obj))

  if(curve.type == "cumulative"){
    sv_list <- lapply(variables, function(x){
	    f <- as.formula(paste("Surv(", time.var, ", ", event.var, ") ~ ", x))
	    survminer::ggsurvplot(survminer::surv_fit(f, data = pop.cutoff.obj[,c(x, time.var, event.var)]),
	               surv.median.line = "hv", pval = TRUE, fun = "event",
	               palette = palette, legend.title = "", title = x,
	               risk.table = TRUE, risk.table.height = 0.25,
	               legend.labs = levels(pop.cutoff.obj[,x]),
	               ggtheme = ggpubr::theme_pubclean())})
  }else{
    sv_list <- lapply(variables, function(x){
	    f <- as.formula(paste("Surv(", time.var, ", ", event.var, ") ~ ", x))
	    survminer::ggsurvplot(survminer::surv_fit(f, data = pop.cutoff.obj[,c(x, time.var, event.var)]),
	               surv.median.line = "hv", pval = TRUE,
	               palette = palette, legend.title = "", title = x,
	               risk.table = TRUE, risk.table.height = 0.25,
	               legend.labs = levels(pop.cutoff.obj[,x]),
	               ggtheme = ggpubr::theme_pubclean())})
  }

  plotaux <- lapply(1:length(sv_list), function(x){
    cowplot::plot_grid(sv_list[[x]]$plot, sv_list[[x]]$table, nrow = 2, rel_heights = c(1,.4))})
  cowplot::plot_grid(plotlist = lapply(plotaux, print))  
}
