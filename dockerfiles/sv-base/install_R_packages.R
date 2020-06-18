#!/usr/bin/env Rscript

# treat warnings as errors, otherwise script can fail silently if a package fails to install
options(warn = 2)

# fixed place to dump downloaded sources
PKG_DUMP="/tmp/R_pkg_download/"

# packages to be installed
to.be.installed = commandArgs(trailingOnly=TRUE)
if ( 0 == length(to.be.installed) ) {
    stop("At least one argument must be supplied when you call me.n", call.=FALSE)
}

# multiple repos, multiple retries when a package is not found
repos <- c("http://lib.stat.cmu.edu/R/CRAN/")

# install to default place, quietly, then leave
install.packages(pkgs = to.be.installed, 
                 repos = repos, 
                 destdir = PKG_DUMP, 
                 clean = TRUE, 
                 quiet = TRUE)
q(save = "no")
