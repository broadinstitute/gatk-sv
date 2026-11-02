#!/usr/bin/env Rscript
# override for debugging purposes
quiet <- as.logical(Sys.getenv("QUIET", unset="TRUE"))
library(parallel)
Ncpus <- detectCores()

# treat warnings as errors, otherwise script can fail silently if a package fails to install
options(warn = 2)

# packages to be installed
to.be.installed <- commandArgs(trailingOnly=TRUE)
if ( 0 == length(to.be.installed) ) {
    stop("At least one argument must be supplied when you call me.n", call.=FALSE)
}

# multiple repos, multiple retries when a package is not found
# NOTE: the old primary mirror (http://lib.stat.cmu.edu/R/CRAN/) now redirects
# plain HTTP to HTTPS / 403s, and with options(warn=2) an unreachable first repo
# aborts the install even when a fallback mirror is listed. Use HTTPS CRAN CDN.
repos <- c("https://cloud.r-project.org", "https://cran.rstudio.com")

# install to default place, quietly, then leave
install.packages(pkgs = to.be.installed, 
                 repos = repos, 
                 clean = TRUE,
                 quiet = quiet,
                 Ncpus = Ncpus)
q(save = "no")
