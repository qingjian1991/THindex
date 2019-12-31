#create a new
#http://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html



#some information about developmets.

library(devtools)

#1) Creating the Framework for your First Package

#devtools::create("THindex")

##2) add dependence
#use_package("dplyr") # Defaults to imports
#> Adding dplyr to Imports
#> Refer to functions with dplyr::fun()
#use_package("dplyr", "Suggests")
use_package("magrittr")
use_package("vegan")
use_package("data.table")
use_package("maftools")
#> Adding dplyr to Suggests
#> Use requireNamespace("dplyr", quietly = TRUE) to test if package is
#>  installed, then use dplyr::fun() to refer to functions.

use_gpl3_license("GPLv3")
use_vignette("THindex")
#5) How do I Document My Functions?
devtools::document()

#4) Edit your code and load your code
devtools::load_all()


devtools::install(build_vignettes = T)



library(THindex)
library(maftools)

laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read.maf(maf = laml.maf)
TCGA.ab.het <- inferHeterogeneityPlus(maf = laml, vafCol = 'i_TumorVAF_WU', index = "diversity")
print(TCGA.ab.het$diveristy)







