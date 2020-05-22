<img src="extra/logo_FlowCT_hex_gihub.png" height="200" align="right" />

# FlowCT: A semi-automated workflow for deconvolution of immunophenotypic data and objective reporting on large datasets 

![](https://img.shields.io/badge/devel%20version-2.1-orange.svg)
[![Last-changedate](https://img.shields.io/badge/last%20change-2020--05--22-green.svg)](https://github.com/jgarces02/FlowCT/commits/devel)

FlowCT is a semi-automated pipeline for flow cytometry data analysis. 
Starting from compensated data obtained with standardized protocols, allows simultaneous analyses of multiple files, automated cell clustering and statistical analysis. It provides results in tabular format that can be exported into other databases for integrated analysis (e.g. clinical trials).

```
options(rsconnect.http = "internal")
Sys.setenv(http_proxy = "http://proxy.unav.es:8080")
Sys.setenv(https_proxy = "http://proxy.unav.es:8080")  
devtools::install_github("jgarces02/FlowCT", auth_token = "21ea9880f944d42755479e54a5b19ddd00fe17f6")
```
