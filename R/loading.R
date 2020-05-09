.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Please, if you use FlowCT package you can cite it as:\n
                            Botta, xxxxxx\n\nAnd for any issue you have, go to our GitHub page: xxx\n")
}

# automatically load required packages when library(FlowCT)
# .onLoad <- function(libname, pkgname) {
#   load_packages <- c("SummarizedExperiment", "flowCore", "ggplot2", "pheatmap", "Rtsne", "uwot", 
#                    "reshape2", "dplyr")
#   lapply(load_packages, suppressPackageStartupMessages(library), character.only = TRUE)
# }
