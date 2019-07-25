.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to my package")
}


.onLoad <- function(libname, pkgname){
  load_packages <- c("readxl", "flowCore", "flowAI", "flowViz", "flowStats", "gridExtra", "ggsci", "matrixStats", "ggplot2", "reshape2", 
                   "ggrepel", "dplyr", "RColorBrewer", "pheatmap", "FlowSOM", "ConsensusClusterPlus", "Rtsne", "uwot", 
                   "premessa", "phytools", "ggtree", "Hmisc", "corrplot", "ggthemes", "ggpubr", "matrixTests", "DataCombine")
  lapply(load_packages, suppressPackageStartupMessages(library), character.only = TRUE)
}
