# 'dotplot.DE
#'
#' It draws dotplot for each cell cluster identified and each marker to facilitate identification of cell populations in each cell cluster.
#' @param fcs.SE A FCS.SE object generated through \code{\link[FlowCT:fcs.SE]{FlowCT::fcs.SE()}}.
#' @param assay.i Name of matrix stored in the \code{FCS.SE} object from which calculate correlation. Default = \code{"normalized"}.
#' @param clusters.named Column name from the \code{initial.fcs.SE} object which contains renamed clusters (through \code{\link[FlowCT:clusters.rename]{FlowCT::clusters.rename()}}).
#' @param markers.to.use Vector with markers to use. Default = \code{"all"}.
#' @param psig.cutoff P-value cutoff. Default = \code{0.05}.
#' @param return.stats Logical indicating if calculated statistics should be returned in a new variable. Default = \code{FALSE}.
#' @param scale.size Numerical value indicating how much scale points. Default = \code{9}.
#' @keywords differential dotplot
#' @keywords cell cluster identification
#' @keywords cell cluster markers
#' @export
#' @import ggplot2
#' @examples
#' \dontrun{
#' dotplot.DE(fcs.SE = fcs_se, markers.to.use = c("CD8", "CD27", "CCR4", "CD45RA", "CD4"), clusters.named = "SOM_named")
#' }

dotplot.DE <- function(fcs.SE, assay.i = "normalized", clusters.named = "SOM_named", markers.to.use = "all", psig.cutoff = 0.05, return.stats = F, scale.size = 9){
  if(markers.to.use == "all") markers.to.use <- rownames(fcs.SE)
  
  dt <- data.frame()
  for(i in unique(fcs.SE[[clusters.named]])){
    aux_se <- fcs.SE[,fcs.SE[[clusters.named]] == i]
    dt <- rbind(dt, data.frame(t(assay(aux_se, i = assay.i)[markers.to.use,]), pop = i))
  }
  dtm <- data.table::melt(dt)
  
  dtm2 <- data.frame()
  for(i in unique(dtm$variable)){
    s1 <- as.data.frame(summary(stats::lm(value ~ 0 + pop, dtm[dtm$variable == i,]))$coefficients)
    s1$pop <- gsub("pop", "", rownames(s1))
    s1$`t value` <- ifelse(s1$`Pr(>|t|)` < psig.cutoff, s1$`t value`, NA) #not consider those non-significant
    s1$mixstats <- s1$Estimate*s1$`t value`
    
    s2 <- stats::aggregate(value ~ ., dtm[dtm$variable == i,], FUN = stats::median) #collapse to median cell values
    s2 <- merge(s2, s1, by = "pop")
    dtm2 <- rbind(dtm2, s2)
    dtm2$pop <- factor(dtm2$pop)
  }
  
  print(ggplot(dtm2, aes_string("variable", y = "pop", size = value, color = "mixstats") + 
    geom_point() + 
    scale_color_continuous(na.value = "gray70", name = "Marker\nimportance") +
    scale_size_area(max_size = scale.size, name = "Median\nfluorescence") +
    labs(x = "", y = "") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  
  if(return.stats) return(dtm2)
}

