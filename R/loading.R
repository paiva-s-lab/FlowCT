.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Please, if you use FlowCT package you can cite it as:\n
                            Botta, xxxxxx\n\nAnd for any issue you have, go to our GitHub page: xxx")
}

# automatically load required packages when library(FlowCT)
# .onLoad <- function(libname, pkgname) {
#   load_packages <- c("flowCore", "flowAI", "flowViz", "flowStats", "gridExtra", "ggsci", "matrixStats", "ggplot2", "reshape2", 
#                      "ggrepel", "dplyr", "RColorBrewer", "pheatmap", "FlowSOM", "ConsensusClusterPlus", "Rtsne", "uwot", 
#                      "premessa", "phytools", "ggtree", "Hmisc", "corrplot", "ggthemes", "ggpubr", "matrixTests", "DataCombine")
#   lapply(load_packages, suppressPackageStartupMessages(library), character.only = TRUE)
# }