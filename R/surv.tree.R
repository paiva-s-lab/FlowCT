
#' Kaplan-Meier from survival tree
#'
#' Plot with a Kaplan-Meier curve with node leafs identified by recursive partitioning.
#' @param rpart.tree A object generated through \code{\link[FlowCT:rpart.tree]{FlowCT::rpart.tree()}}.
#' @param tree Specify the tree name, i.e. "original" or "pruned".
#' @param time.var Survival time variable.
#' @param event.var Variable with event censoring.
#' @param return.data Metadata with node leafs distribution should be returned. Default = \code{FALSE}.
#' @param curve.type Plotting curve methodology, i.e. "survival" (default) or "cumulative" (for progressions).
#' @keywords recursive partioning
#' @keywords CART
#' @keywords classification regression tree
#' @keywords survival cutoffs Kaplan-Meier curves
#' @export rpart.tree
#' @examples
#' \dontrun{
#' surv.tree(rpart.tree = tr, tree = "original", time.var = "PFS", event.var = "PFS_c")
#' }

surv.tree <- function(rpart.tree, tree, time.var, event.var, return.data = F, curve.type = "survival"){
  if (!requireNamespace(c("survminer", "partykit", "cowplot"), quietly = TRUE)) stop("Packages \"survminer\", \"cowplot\" and \"partykit\" needed for this function to work. Please install them.", call. = FALSE)

  md <- rpart.tree$metadata
  tr <- rpart.tree[[grep(tree, names(rpart.tree))]]
  tr$splits <- round(tr$splits, 2)

  if(class(md[,event.var]) != "numeric") md[,event.var] <- as.numeric(as.factor(md[,event.var]))
  if(class(md[,time.var]) != "numeric") md[,time.var] <- as.numeric(md[,time.var])

  md$group <- tr$where
  
  nodelabs <- partykit:::.list.rules.party(partykit::as.party.rpart(tr, data = TRUE))
  nodelabs <- sapply(nodelabs, function(i){ #round cutoffs
    aux <- strsplit(i, " ")[[1]]
    suppressWarnings(paste(ifelse(!is.na(as.numeric(aux)), round(as.numeric(aux), 2), aux), collapse = " "))
  })

  md$group_cod <- factor(md$group)
  levels(md$group_cod) <- nodelabs

  f <- as.formula(paste0("Surv(", time.var, ", ", event.var, ") ~ group"))
  if(curve.type == "cumulative"){
    g1 <- survminer::ggsurvplot(survminer::surv_fit(f, data = md), pval = T, surv.median.line = "hv",
               ggtheme = theme_light(), legend.labs = nodelabs, fun = "event",
               risk.table = T, risk.table.y.text = F) + guides(colour = guide_legend(ncol = 1))
  }else{
    g1 <- survminer::ggsurvplot(survminer::surv_fit(f, data = md), pval = T, surv.median.line = "hv",
               ggtheme = theme_light(), legend.labs = nodelabs,
               risk.table = T, risk.table.y.text = F) + guides(colour = guide_legend(ncol = 1))
  }

  print(g2 <- cowplot::plot_grid(g1$plot, g1$table, nrow = 2, rel_heights = c(1,.4)))

  if(return.data) return(list(metadata = md, plot = g2))
}
