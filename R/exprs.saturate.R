#' exprs.saturate
#' @export

exprs.saturate <- function(expression_df){
  rng <- matrixStats::colQuantiles(expression_df, probs = c(0.01, 0.99))
  expr01 <- t((t(expression_df) - rng[, 1]) / (rng[, 2] - rng[, 1]))
  expr01[expr01 < 0] <- 0
  expr01[expr01 > 1] <- 1
  
  return(expr01)
}
