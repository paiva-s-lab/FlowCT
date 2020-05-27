#' div.colors
#' @export

div.colors <- function(n){
  require(RColorBrewer)

  if (n < 74) { # maximum number of unique colors combining all palettes
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual", ]
    qual_col_pals = qual_col_pals[c("Paired", "Set1", "Dark2", "Accent", "Set2", "Set3", "Pastel1", "Pastel2"),]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col = unique(col_vector)[1:n]
  }else{ # up to 295 colors without repetitive colors
    qual_col_pals = sample(brewer.pal.info)
    col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
    set.seed(333); col <- col_vector[unique(sample(n))]
  }
  return(col)
}
