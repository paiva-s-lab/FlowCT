.onAttach <- function(libname, pkgname) {
  require(SingleCellExperiment)
  packageStartupMessage("\n --- --- --- \n
    Please, if you use FlowCT package you can cite it as:\n
    Cirino Botta, Catarina Da Silva Maia, Juan-José Garcés et al.
    FlowCT for the analysis of large immunophenotypic datasets and biomarker discovery in cancer immunology.
    Blood Advances 2021.
    \n --- --- --- \n")
}
