### ---------
#' Convert table to ASCII format
#' @export

asciify <- function(df, pad = 1, ...) {
    ## error checking
    stopifnot(is.data.frame(df))
    ## internal functions
    SepLine <- function(n, pad = 1) {
        tmp <- lapply(n, function(x, pad) paste(rep("-", x + (2* pad)),
                                                collapse = ""), pad = pad)
        paste0("+", paste(tmp, collapse = "+"), "+")
    }
    Row <- function(x, n, pad = 1) {
        foo <- function(i, x, n) {
            fmt <- paste0("%", n[i], "s")
            sprintf(fmt, as.character(x[i]))
        }
        rowc <- sapply(seq_along(x), foo, x = x, n = n)
        paste0("|", paste(paste0(rep(" ", pad), rowc, rep(" ", pad)),
                          collapse = "|"), "|")
    }
    
    ## convert everything to characters
    df <- as.matrix(df)
    ## nchar in data
    mdf <- apply(df, 2, function(x) max(nchar(x)))
    ## nchar in names
    cnames <- nchar(colnames(df))
    ## max nchar of name+data per elements
    M <- pmax(mdf, cnames)

    if(nrow(df) < 20){
      ## write the header
      sep <- SepLine(M, pad = pad)
      writeLines(sep)
      writeLines(Row(colnames(df), M, pad = pad))
      writeLines(sep)
      ## write the rows
      for(i in seq_len(nrow(df))) {
          ## write a row
          writeLines(Row(df[i,], M, pad = pad))
          ## write separator
          writeLines(sep)
      }
    }else{
      sep <- SepLine(M, pad = pad)
      writeLines(Row(colnames(df), M, pad = pad))
      writeLines(sep)
      for(i in seq_len(nrow(df))) {
          writeLines(Row(df[i,], M, pad = pad))
      }
    }
    invisible(df)
}

### ---------
#' Create a palette of colors
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

### ---------
#' Create a palette of colors suitable for doing a heatmap
#' @export

col.annot.pheatmap <- function(metadata, colors = NULL){
  annot_col <- list()
  for(i in colnames(metadata)){
    if(is.numeric(metadata[[i]])) next
    aux <- as.vector(unlist(unique(metadata[i])))
    if(is.null(colors)) colors2 <- div.colors(length(aux)) else colors2 <- colors
    names(colors2) <- aux
    annot_col[[i]] <- colors2
  }
  return(lapply(annot_col, function(x) x[!is.na(names(x))])) # delete NAs for list with less elements specified in "colors"
}

### ---------
#' Transform a \code{fcs.SCE} object into a \code{flowset} object
#' @export

as.flowSet.SE <- function(fcs.SCE, assay.i){
  df_flowset <- lapply(unique(fcs.SCE$filename), function(x) as_flowFrame.me(t(SummarizedExperiment::assay(fcs.SCE[,fcs.SCE$filename == x], i = assay.i))))
  names(df_flowset) <- unique(fcs.SCE$filename)
  return(as(df_flowset, "flowSet"))
}

# modified from https://github.com/ParkerICI/premessa/blob/master/R/fcs_io.R
as_flowFrame.me <- function(exprs.m) {
    flow.frame <- flowCore::flowFrame(exprs.m)
	params <- flowCore::parameters(flow.frame)
    pdata <- flowCore::pData(params)

    for (i in 1:ncol(flow.frame)) {
        s <- paste("$P",i,"S",sep="")
        n <- paste("$P",i,"N",sep="")
        r <- paste("$P",i,"R",sep="")
        b <- paste("$P",i,"B",sep="")
        e <-  paste("$P",i,"E",sep="")

        keyval <- list()

        keyval[[s]] <- colnames(exprs.m)
        keyval[[n]] <- colnames(exprs.m)[i]
        keyval[[r]] <- ceiling(max(exprs.m[,i], na.rm = TRUE))

        keyval[[b]] <- 32
        keyval[[e]] <- "0,0"
        flowCore::keyword(flow.frame) <- keyval

        pdata[i,"minRange"] <- min(exprs.m[,i], na.rm = TRUE)
        pdata[i,"maxRange"] <- max(exprs.m[,i], na.rm = TRUE)
    }

    flowCore::pData(params) <- pdata
    flowCore::parameters(flow.frame) <- params

    return(flow.frame)
}
