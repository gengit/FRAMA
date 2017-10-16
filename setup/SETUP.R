#!/usr/bin/env Rscript

# --------------------------------------------------------------------
# Install required packages from CRAN and Bioconductor
#
# Martin Bens | bensmartin@gmail.com
# 2017-08-27
# --------------------------------------------------------------------

personal.lib = Sys.getenv("R_LIBS_USER")
if (!file.exists(personal.lib)) { dir.create(personal.lib, rec = T) }

# CRAN packages
mypackages = c("plyr", "ggplot2", "reshape", "gridExtra", "grid")
installme = mypackages[!(mypackages %in% installed.packages()[,"Package"])]
if (length(installme) > 0) {
    paste("The following packages need to be installed:", paste(installme, collapse = ", "))
    install.packages(installme, repos =  "http://cran.us.r-project.org", lib = personal.lib)
}

# Bioconductor packages
source("http://bioconductor.org/biocLite.R")
mypackages = c("GO.db", "org.Hs.eg.db", "KEGG.db")
installme = mypackages[!(mypackages %in% installed.packages()[,"Package"])]
if(length(installme) > 0) {
    paste("The following packages need to be installed:", paste(installme, collapse = ", "))
    paste("")
    biocLite(installme, suppressAutoUpdate = T, suppressUpdates = T, lib = personal.lib)
}
