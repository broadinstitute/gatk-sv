#!/usr/bin/env Rscript
# override for debugging purposes
quiet <- as.logical(Sys.getenv("QUIET", unset="TRUE"))

# treat warnings as errors, otherwise script can fail silently if a package fails to install
options(warn = 2)

# packages to be installed
to.be.installed <- commandArgs(trailingOnly=TRUE)
if ( 0 == length(to.be.installed) ) {
    stop("At least one argument must be supplied when you call me.n", call.=FALSE)
}

# multiple repos, multiple retries when one of them is un-reachable
repos <- c("https://cran.rstudio.com",
           "http://cran.mtu.edu",
           "http://lib.stat.cmu.edu/R/CRAN/")

# install to default place, quietly, then leave
if (!requireNamespace("BiocManager")){
    install.packages(pkgs = "BiocManager", 
                     repos = repos, 
                     clean = TRUE,
                     quiet = quiet)
}
BiocManager::install(pkgs = to.be.installed, 
                     ask = FALSE, 
                     clean = TRUE,
                     quiet = quiet)
q(save = "no")
