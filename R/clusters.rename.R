#' Rename clusters from numeric to cell name
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
#' fcs$SOM_named <- clusters.rename(fcs$SOM, cluster = replacedata$original_cluster, 
#' 		name = replacedata$new_cluster)
#' 
#' # option 2
#' fcs$SOM_named <- clusters.rename(fcs$SOM, cluster = 1:6, 
#' 		name = c("debris", "lymphocytes", "monocytes", "lymphocytes", "eosinophils", 
#' 		"eosinophils"))
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
