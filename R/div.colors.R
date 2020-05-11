#' div.colors
div.colors <- function(n, set.seed = 333){
  require(RColorBrewer)
  
  if(n < 74){ #all posible colors in qual type palettes
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col = unique(col_vector)
    set.seed(set.seed); col <- sample(col, n, replace = T)
  }else{
    qual_col_pals = sample(brewer.pal.info)
    col_vector = sample(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
    col = sample(col_vector, 333)
    set.seed(set.seed); col <- sample(col, n, replace = T)
  }
  return(col)
}
