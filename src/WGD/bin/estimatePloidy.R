#!/usr/bin/env Rscript

# Script to analyze a matrix of coverage (binCov) values per contig per sample.
# Performs the following:
# 1) copy number estimation per contig per sample (with visualization options)
# 2) sex assignment (with visualization options)
# 3) [optional] PCA & k-means based sample batching (with visualization options)

#Input: a binCov matrix of any resolution (very coarse resolution recommended, e.g. 1Mb)

####################################
#####Set parameters & load libraries
####################################
options(scipen=1000, stringsAsFactors=F, check.names=F)

############################################
#####Helper function to load coverage matrix
############################################
readMatrix <- function(INFILE){
  dat <- read.table(INFILE, comment.char="", header=T, check.names=F)
  colnames(dat)[1] <- "Chr"
  dat[, -1] <- t(apply(dat[, -1], 1, as.numeric))
  return(dat)
}

#############################################################
#####Helper function to normalize contigs for a single sample
#############################################################
normalizeContigsPerSample <- function(mat, exclude=c("X", "Y"), ploidy=2){
  #Convert vals to numeric
  mat[, -c(1:3)] <- apply(mat[, -c(1:3)], 2, as.numeric)

  #Iterate over all samples and scale each sample by median and center/scale to expected ploidy
  mat[, -c(1:3)] <- sapply(4:ncol(mat), function(i){
    #Compute median of values excluding current value & any other specified
    sampVals <- as.vector(mat[which(!(mat[, 1] %in% exclude)), i])
    excl.median <- median(sampVals[which(sampVals>0)], na.rm=T)

    #Normalize values by excl.median
    newVals <- mat[, i]/excl.median

    #Scale to expected ploidy
    newVals <- ploidy*newVals

    #Return cleaned values
    return(newVals)
  })

  #Return normalized matrix
  return(mat)
}

###########################################################################
#####Helper function to remove rows where >X% of samples have <=Y% coverage
###########################################################################
filterZeroBins <- function(mat, exclude=c("X", "Y"), minSamp=0.8, minCov=0.2){
  #Convert vals to numeric
  mat[, -c(1:3)] <- apply(mat[, -c(1:3)], 2, as.numeric)

  #Find bins where >minSamp% of samples have coverage>minCov
  #Dev note, Jan 2020: default for minSamp boosted from 0.05 to 0.8, 
  # and minCov boosted from 0.05 to 0.2, to account for switch 
  # from binCov to GATK collectReadCounts
  nZeros <- apply(mat[, -c(1:3)], 1, function(vals){
    return(length(which(vals<=minCov)))
  })
  fracZeros <- nZeros/(ncol(mat)-3)
  keep.bins <- which(fracZeros<=minSamp | mat[, 1] %in% exclude)

  #Return subsetted matrix
  return(mat[keep.bins, ])
}

##################################################
#####Helper function to run PCA on a binCov matrix
##################################################
binCovPCA <- function(dat, exclude=c("X", "Y"), topPCs=10){
  #Runs PCA
  PCA <- prcomp(dat[which(!(dat[, 1] %in% exclude)), -c(1:3)], center=T, scale=T)

  #Restricts to the loadings for the top N PCs
  PCs <- as.data.frame(PCA$rotation[, 1:topPCs])

  #Formats & returns PC data frame
  out.df <- data.frame("sample_id"=names(dat[, -c(1:3)]))
  out.df <- cbind(out.df, PCs)
  rownames(out.df) <- 1:nrow(out.df)
  return(list("full"=PCA, "top"=out.df))
}

############################################################
#####Helper function to compute median per contig per sample
############################################################
medianPerContigPerSample <- function(dat){
  #Get list of unique contigs
  contigs <- unique(dat[, 1])

  #Iterate over all contigs and compute median bin value per sample
  allMedians <- sapply(contigs, function(contig){
    #Calculate median value per sample
    chrVals <- apply(dat[, -c(1:3)], 2, function(vals){
      #Subset vals to contig
      vals <- vals[which(dat[, 1]==contig)]
      #Return median
      return(median(vals, na.rm=T))
    })
    #Return vector of medians
    return(as.vector(chrVals))
  })

  #Compose output data frame
  out.df <- data.frame("sample_id"=names(dat[, -c(1:3)]))
  out.df <- cbind(out.df, allMedians)
}

########################################################
#####Helper function to PCA-cluster samples into batches
########################################################
clusterDosageProfiles <- function(PCs, n.min=60, n.max=150, n.ideal=100, wait.max=10, plot=T){
  #Get desired k & target number of samples per cluster
  k <- round(nrow(PCs)/n.ideal)
  k.colors <- rainbow(k)
  n.target <- round(nrow(PCs)/k)

  #Reformat PC data (assumes first column is ID)
  rownames(PCs) <- PCs[, 1]
  PCs <- as.data.frame(PCs[, -1])

  #Initial k-means clustering
  clust.init <- kmeans(PCs, centers=k)
  clust <- clust.init
  centers.init <- clust$centers
  centers.prev <- centers.init
  centers <- centers.init

  #Initialize PC1 vs PC2 plot
  if(plot==T){
    layout(matrix(c(1, 2, 5, 3, 4, 5), nrow=2, byrow=T), widths=c(3, 3, 1))
    plot(PCs[, 1], PCs[, 2], pch=19, cex=0.2, 
         main="Center Progression", xlab="PC1", ylab="PC2")
    points(centers.init[, 1], centers.init[, 2], pch=21, bg=k.colors)
  }

  #Adjust cluster centers
  #Loop until no clusters are > n.max or are < n.min, 
  # or cluster centers stop moving after at least wait.max tries
  j=0
  while(any(clust$size<n.min | clust$size>n.max) & !(j>=wait.max & all(centers-centers.prev==0))){
    #Calculate distance between each sample and cluster centers
    samp.dist <- as.data.frame(t(apply(PCs, 1, function(vals){
      #Iterate over clusters
      apply(centers, 1, function(clusterCoords){
        dist(rbind(vals, clusterCoords))
      })
    })))

    #Iterate over all clusters, add/subtract samples up to n.min/n.max, and recompute cluster center
    centers <- as.data.frame(t(sapply(1:nrow(centers), function(i){
      #Check if cluster is too small
      if(clust$size[i]<n.min){
        #Get closest n.min samples
        nuc <- head(order(samp.dist[, i]), n.min)
        #Compute new cluster center coordinates & assign to centers df
        new.cent <- apply(PCs[nuc, ], 2, mean)
      }else{
        #Check if cluster is too large
        if(clust$size[i]>n.max){
          #Get closest n.max samples
          nuc <- head(order(samp.dist[, i]), n.max)
          #Compute new cluster center coordinates & assign to centers df
          new.cent <- apply(PCs[nuc, ], 2, mean)
        }else{
          #Otherwise, recompute center according to best n.target assigned samples
          nuc <- head(order(samp.dist[, i]), n.target)
          #Compute new cluster center coordinates & assign to centers df
          new.cent <- apply(PCs[nuc, ], 2, mean)
        }
      }
      return(new.cent)
    })))

    #Rerun k-means with new centers
    clust <- kmeans(PCs, centers=centers)
    centers <- clust$centers

    #Update plot, if optioned
    if(plot==T){
      segments(x0=centers.prev[, 1], x1=centers[, 1], 
               y0=centers.prev[, 2], y1=centers[, 2], 
               col=k.colors)
      points(centers[, 1], centers[, 2], pch=21, bg=k.colors)
    }
    centers.prev <- centers

    #Increment counter
    j <- j+1
  }

  #Update which clusters are small/large
  large.clusters <- which(clust$size>n.max)
  small.clusters <- which(clust$size<n.min)

  #Calculate distance between each sample and cluster centers
  samp.dist <- as.data.frame(t(apply(PCs, 1, function(vals){
    #Iterate over clusters
    apply(centers, 1, function(clusterCoords){
      dist(rbind(vals, clusterCoords))
    })
  })))

  ####Procedure to make obvious switches between large and small clusters
  if(length(large.clusters)>0 & length(small.clusters)>0){
    #Get list of samples with poor fits to large clusters
    poor.fits <- unique(as.vector(sapply(large.clusters, function(i){
      #Get all members in cluster
      members <- which(clust$cluster==i)
      #Order members by goodness of fit
      members <- members[as.vector(order(PCs[as.vector(members), i]))]
      #Return tail of members that are poorest fits to large cluster
      as.vector(tail(members, clust$size[i]-n.max))
    })))

    #Get samples that are good fits for one (or more) small clusters
    good.fits <- unique(as.vector(sapply(small.clusters, function(i){
      return(as.vector(head(order(samp.dist[, i]), n.min)))
    })))

    #Reassign poor fits to good fits, if any
    if(any(good.fits %in% poor.fits)){
      for(i in good.fits[which(good.fits %in% poor.fits)]){
        #Get old assignment
        old <- clust$cluster[i]
        #Get best fit among small clusters
        new <- small.clusters[which(samp.dist[i, small.clusters]==min(samp.dist[i, small.clusters]))]
        #Modify cluster assignment
        clust$cluster[i] <- new
      }
    }

    #Recompute cluster size & cluster centers
    clust$size <- as.vector(table(clust$cluster))
    clust$centers <- t(sapply(1:nrow(clust$centers), function(i){
      apply(PCs[which(clust$cluster==i), ], 2, mean)
    }))
    centers <- clust$centers

    #Update list of clusters < n.min & clusters > n.max
    large.clusters <- which(clust$size>n.max)
    small.clusters <- which(clust$size<n.min)

    #Update plot, if optioned
    if(plot==T){
      segments(x0=centers.prev[, 1], x1=centers[, 1], 
               y0=centers.prev[, 2], y1=centers[, 2], 
               col=k.colors)
      points(centers[, 1], centers[, 2], pch=21, bg=k.colors)
    }
    centers.prev <- centers
  }

  ####Procedure to reassign samples such that all clusters are > n.min
  if(length(small.clusters)>0){
    #Calculate distance between each sample and cluster centers
    samp.dist <- as.data.frame(t(apply(PCs, 1, function(vals){
      #Iterate over clusters
      apply(centers, 1, function(clusterCoords){
        dist(rbind(vals, clusterCoords))
      })
    })))

    #Get list of best n.min samples for each cluster with ≥ n.min
    # and all currently assigned samples for each cluster < n.min
    #These samples cannot be reassigned in the subsequent step
    no.swap <- unique(as.vector(sapply(which(clust$size>=n.min), function(i){
      #Get samples assigned to cluster
      members <- which(clust$cluster==i)
      #Order members by distance to cluster center & return top n.min
      return(as.vector(head(members[order(samp.dist[members, i])], n.min)))
    })))
    no.swap <- unique(c(no.swap, which(clust$cluster %in% small.clusters)))

    #Rank samples available for reassignment in terms of best fit for each cluster < n.min
    reassign.ranks <- apply(as.data.frame(samp.dist[-no.swap, small.clusters]), 2, rank, ties.method="random")
    reassign.ranks <- reassign.ranks[order(apply(reassign.ranks, 1, min)), ]

    #Iterate over ranked reassignable samples
    #For each sample, check if the preferred fit is still < n.min
    #If so, reassign it to that cluster. If not, skip it
    #After each sample, recompute cluster sizes
    if(nrow(reassign.ranks)>0){
      l <- 1
      while(any(clust$size<n.min) & k<nrow(reassign.ranks)){
        #Get original sample index
        idx <- which(rownames(PCs)==rownames(reassign.ranks)[l])
        #Get preferred assignment for sample
        pref <- sample(small.clusters[which(reassign.ranks[l, ]==min(reassign.ranks[l, ]))], 1)
        #Check if that cluster still needs samples
        if(clust$size[pref]<n.min){
          #If yes, reassign
          clust$cluster[idx] <- pref
        }
        #Recompute cluster sizes
        clust$size <- as.vector(table(clust$cluster))
        #Update counter
        l <- l+1
      }
    }

    #Recompute cluster size & cluster centers
    clust$size <- as.vector(table(clust$cluster))
    clust$centers <- t(sapply(1:nrow(clust$centers), function(i){
      apply(PCs[which(clust$cluster==i), ], 2, mean)
    }))
    centers <- clust$centers

    #Update plot, if optioned
    if(plot==T){
      segments(x0=centers.prev[, 1], x1=centers[, 1], 
               y0=centers.prev[, 2], y1=centers[, 2], 
               col=k.colors)
      points(centers[, 1], centers[, 2], pch=21, bg=k.colors)
    }
    centers.prev <- centers
  }

  #Recompute cluster size & cluster centers
  clust$size <- as.vector(table(clust$cluster))
  clust$centers <- t(sapply(1:nrow(clust$centers), function(i){
    apply(PCs[which(clust$cluster==i), ], 2, mean)
  }))
  centers <- clust$centers

  #Recompute sample distance to cluster centers
  samp.dist <- as.data.frame(t(apply(PCs, 1, function(vals){
    #Iterate over clusters
    apply(centers, 1, function(clusterCoords){
      dist(rbind(vals, clusterCoords))
    })
  })))

  #Update list of clusters > n.max
  large.clusters <- which(clust$size>n.max)
  small.clusters <- which(clust$size<n.min)

  ####Procedure to reassign samples such that all clusters are < n.max
  if(length(large.clusters>0)){
    #Get list of worst samples for each cluster with ≥ n.min
    # such that each cluster retains exactly n.max samples
    #These samples will be reassigned in the subsequent step
    swappable <- unique(as.vector(unlist(sapply(large.clusters, function(i){
      #Get samples assigned to cluster
      members <- which(clust$cluster==i)
      #Order members by distance to cluster center & return worst samples after skipping n.max
      return(as.vector(tail(members[order(samp.dist[members, i])], clust$size[i]-n.max)))
    }))))

    #Rank samples available for reassignment in terms of best fit for each cluster < n.min
    if(length(swappable)>0){
      reassign.ranks <- apply(as.data.frame(samp.dist[swappable, -large.clusters]), 2, rank, ties.method="random")
      reassign.ranks <- reassign.ranks[order(apply(reassign.ranks, 1, min)), ]
    }

    #Iterate over ranked reassignable samples
    #For each sample, reassign to top ranked preferred fit
    # as long as reassignment does bring preferred fit to > n.max
    #After each sample, recompute cluster sizes
    if(nrow(reassign.ranks)>0){
      for(l in 1:nrow(reassign.ranks)){
        #Get original sample index
        idx <- which(rownames(PCs)==rownames(reassign.ranks)[l])
        #Rank possible sample assignments by distance
        pref <- (1:nrow(centers))[-large.clusters][order(reassign.ranks[l, ])]
        #Assign to top preferred cluster, provided that cluster isn't > n.max
        new <- head(pref[pref %in% which(clust$size<n.max)], 1)
        clust$cluster[idx] <- new
        #Recompute cluster sizes
        clust$size <- as.vector(table(clust$cluster))
      }
    }

    #Recompute cluster size & cluster centers
    clust$size <- as.vector(table(clust$cluster))
    clust$centers <- t(sapply(1:nrow(clust$centers), function(i){
      apply(PCs[which(clust$cluster==i), ], 2, mean)
    }))
    centers <- clust$centers

    #Update plot, if optioned
    if(plot==T){
      segments(x0=centers.prev[, 1], x1=centers[, 1], 
               y0=centers.prev[, 2], y1=centers[, 2], 
               col=k.colors)
      points(centers[, 1], centers[, 2], pch=21, bg=k.colors)
    }
    centers.prev <- centers
  }

  #Recompute cluster size & cluster centers
  clust$size <- as.vector(table(clust$cluster))
  clust$centers <- t(sapply(1:nrow(clust$centers), function(i){
    apply(PCs[which(clust$cluster==i), ], 2, mean)
  }))
  centers <- clust$centers

  #Recompute sample distance to cluster centers
  samp.dist <- as.data.frame(t(apply(PCs, 1, function(vals){
    #Iterate over clusters
    apply(centers, 1, function(clusterCoords){
      dist(rbind(vals, clusterCoords))
    })
  })))

  ####Final optimization
  #Iterate wait.max times
  j <- 0
  while(j<wait.max){
    #Rank matches for each sample vs all clusters
    final.ranks <- as.data.frame(apply(samp.dist, 2, rank, ties.method="random"))
    final.ranks <- final.ranks[order(apply(final.ranks, 1, min)), ]
    final.ranks <- as.data.frame(t(apply(final.ranks, 1, rank, ties.method="random")))
    #Subset to samples that aren't already assigned to their top preference
    samples.to.switch <- sapply(1:nrow(final.ranks), function(i){
      #Get original sample index
      idx <- which(rownames(PCs)==rownames(final.ranks)[i])
      #Get current & preferred assignment for sample
      current <- clust$cluster[idx]
      pref <- which(final.ranks[i, ]==min(final.ranks[i, ]))
      #Return sample if pref != current
      if(pref!=current){
        return(i)
      }else{
        return(NA)
      }
    })
    samples.to.switch <- samples.to.switch[which(!is.na(samples.to.switch))]
    final.ranks <- final.ranks[samples.to.switch, ]
    #Iterate over all samples that aren't already assigned to
    # their top preference. If sample fits another cluster better, 
    # and wouldn't bring the current cluster < n.min, 
    # or the new cluster > n.max, reassign that sample
    if(nrow(final.ranks)>0){
      for(l in 1:nrow(final.ranks)){
        #Get original sample index
        idx <- which(rownames(PCs)==rownames(final.ranks)[l])
        #Get current assignment for sample
        current <- clust$cluster[idx]
        #Try to assign to any cluster with higher priority than current
        if(final.ranks[l, current]>1){
          #Only try if current cluster > n.min
          if(clust$size[current] > n.min){
            #Get ordered list of all better options
            pref <- order(final.ranks[l, ])[which(final.ranks[l, ]<final.ranks[l, current])]
            for(t in pref){
              #Reassign if preferred cluster < n.max
              if(clust$size[t] < n.max){
                clust$cluster[idx] <- t
                break
              }
            }
          }
        }
        #Recompute cluster sizes after each sample
        clust$size <- as.vector(table(clust$cluster))
      }
    }

    #Recompute cluster size & cluster centers
    clust$size <- as.vector(table(clust$cluster))
    clust$centers <- t(sapply(1:nrow(clust$centers), function(i){
      apply(PCs[which(clust$cluster==i), ], 2, mean)
    }))
    centers <- clust$centers

    #Recompute sample distance to cluster centers
    samp.dist <- as.data.frame(t(apply(PCs, 1, function(vals){
      #Iterate over clusters
      apply(centers, 1, function(clusterCoords){
        dist(rbind(vals, clusterCoords))
      })
    })))

    #Update counter
    j <- j+1
  }

  #Recompute cluster size & cluster centers
  clust$size <- as.vector(table(clust$cluster))
  clust$centers <- t(sapply(1:nrow(clust$centers), function(i){
    apply(PCs[which(clust$cluster==i), ], 2, mean)
  }))
  centers <- clust$centers

  #Update plot, if optioned
  if(plot==T){
    segments(x0=centers.prev[, 1], x1=centers[, 1], 
             y0=centers.prev[, 2], y1=centers[, 2], 
             col=k.colors)
    points(centers[, 1], centers[, 2], pch=23, cex=2, bg=k.colors)
  }
  centers.prev <- centers

  #Master PCA plot at end of first three PCs colored by cluster
  if(plot==T){
    plot(PCs$PC1, PCs$PC2, col=k.colors[clust$cluster], pch=19, cex=0.75, 
         main="PC1 vs. PC2", xlab="PC1", ylab="PC2")
    points(centers[, 1], centers[, 2], pch=23, bg=k.colors, cex=2)
    plot(PCs$PC1, PCs$PC3, col=k.colors[clust$cluster], pch=19, cex=0.75, 
         main="PC1 vs. PC3", xlab="PC1", ylab="PC3")
    points(centers[, 1], centers[, 3], pch=23, bg=k.colors, cex=2)
    plot(PCs$PC2, PCs$PC3, col=k.colors[clust$cluster], pch=19, cex=0.75, 
         main="PC2 vs. PC3", xlab="PC2", ylab="PC3")
    points(centers[, 2], centers[, 3], pch=23, bg=k.colors, cex=2)
    #Legend
    par(bty="n", mar=c(4, 0.2, 4, 0.2))
    plot(x=rep(1, times=length(clust$size)), y=1:length(clust$size), 
         xlim=c(0.5, 4), ylim=c(0, length(clust$size)+1), 
         col=k.colors, pch=19, cex=2, 
         xaxt="n", yaxt="n", xlab="", ylab="")
    text(x=rep(1, times=length(clust$size)), y=1:length(clust$size), 
         pos=4, labels=paste("N=", clust$size, sep=""))
  }

  #Return cluster assignments & centers
  out <- list("batch"=clust$cluster, "cluster.centers"=clust$centers)
  return(out)
}

#########################################################################
#####Helper function to normalize contigs for an entire matrix of samples
#########################################################################
normalizeContigsPerMatrix <- function(dat, exclude=NA, scale.exclude=NA, 
                                      genome.ploidy=2, contig.ploidy){
  #Iterate over samples & normalize
  suppressWarnings(if(is.na(exclude)){
    dat[, -1] <- t(apply(dat[, -1], 1, normalizeContigsPerSample, genome.ploidy))
  }else{
    dat[, -1] <- t(apply(dat[, -1], 1, normalizeContigsPerSample, exclude=exclude-1, genome.ploidy))
  })

  #Scale mad to mad of first 12 chromosomes
  mad.others <- mad(unlist(dat[, 2:13]), na.rm=T)

  #Iterate over contigs (minus scale.exclude) and scale
  scaledVals <- sapply(setdiff(2:ncol(dat), scale.exclude), function(i){
    #Calculate & apply adjustments
    median.adjust <- median(dat[, i], na.rm=T)-contig.ploidy[i-1]
    newvals <- dat[, i]-median.adjust
    return(newvals)
  })
  suppressWarnings(if(is.na(scale.exclude)){
    dat[, -1] <- scaledVals
  }else{
    dat[, -c(1, scale.exclude)] <- scaledVals
  })

  #Round up values that were normalized below zero
  dat[, -1] <- apply(dat[, -1], 2, function(vals){
    vals[which(vals<0 & !is.na(vals))] <- 0
    return(vals)
  })

  #Return transformed data
  return(dat)
}

#######################################################################
#####Helper function to test each contig per sample for evidence of CNA
#######################################################################
testCNs <- function(dat, FDR=T){
  dat.out <- dat

  #Iterate over contigs
  left.p <- apply((apply(dat[, -1], 2, scale, center=T, scale=T)), 2, pnorm)
  right.p <- apply((apply(dat[, -1], 2, scale, center=T, scale=T)), 2, pnorm, lower.tail=F)

  #Choose minimum p-value between left and right tails
  dat.out[, -1] <- t(sapply(1:nrow(left.p), function(row){
    pvals <- sapply(1:ncol(left.p), function(col){
      return(min(left.p[row, col], right.p[row, col], na.rm=T))
    })
    return(pvals)
  }))

  #FDR correct all p-values (if optioned)
  if(FDR==T){
    dat.out[, -1] <- apply(dat.out[, -1], 2, p.adjust, method="fdr")
  }

  #Return
  return(dat.out)
}

###########################################
#####Helper function to assign & plot sexes
###########################################
assignSex <- function(dat, sexChr=24:25, 
                      sexAssign.df, #four-column df: cn(X), cn(Y), label, color
                      mosaicThresh="Bonferroni", 
                      plot=T, axLim=3, highlightSample=""){
  #Exclude incomplete entries
  sample.exclude <- unlist(sapply(1:nrow(dat), function(i){
    if(any(is.na(dat[i, sexChr]))){
      return(i)
    }
  }))
  if(length(sample.exclude)>0){
    dat.mod <- dat[-sample.exclude, c(1, sexChr)]
    if(nrow(dat)!=nrow(dat.mod)){
      warning(paste(nrow(dat)-nrow(dat.mod), 
                    " samples missing sex chromosome coverage information and were excluded from sex assignments.", 
                    sep=""))
    }
  }else{
    dat.mod <- dat[, c(1, sexChr)]
  }

  #Round X & Y CN-types to nearest whole integer
  dat.mod.rounded <- dat.mod
  dat.mod.rounded[, -1] <- apply(dat.mod[, -1], 2, round, digits=0)

  #Create output data frame with sex assignments
  sexes <- as.data.frame(t(unlist(sapply(unique(dat[, 1]), function(sample_id){
    #Iterate over all IDs
    if(sample_id %in% dat.mod[, 1]){
      #Get rounded CNs
      CN.X <- round(dat.mod.rounded[which(dat.mod.rounded$sample_id==sample_id), 2], 0)
      CN.Y <- round(dat.mod.rounded[which(dat.mod.rounded$sample_id==sample_id), 3], 0)

      #Assign to existing sex copy-type, or "OTHER"
      if(length(which(sexAssign.df$CN.X==CN.X & sexAssign.df$CN.Y==CN.Y))==1){
        return(as.vector(c(sample_id, CN.X, CN.Y, 
                           sexAssign.df[which(sexAssign.df$CN.X==CN.X & sexAssign.df$CN.Y==CN.Y), 3])))
      }else{
        return(as.vector(c(sample_id, CN.X, CN.Y, "OTHER")))
      }
    }else{
      return(as.vector(c(sample_id, 
                         dat[which(dat$sample_id==sample_id), sexChr[1]], 
                         dat[which(dat$sample_id==sample_id), sexChr[2]], 
                         NA)))
    }
  }))))
  colnames(sexes) <- c("sample_id", "chrX.CN", "chrY.CN", "Assignment")
  rownames(sexes) <- 1:nrow(sexes)

  #Gather sd of X and Y assignments from males
  if(length(sample.exclude)>0){
    sd.X <- sd(dat.mod[which(sexes$Assignment[-sample.exclude]=="MALE"), 2], na.rm=T)
    sd.Y <- sd(dat.mod[which(sexes$Assignment[-sample.exclude]=="MALE"), 3], na.rm=T)
  }else{
    sd.X <- sd(dat.mod[which(sexes$Assignment=="MALE"), 2], na.rm=T)
    sd.Y <- sd(dat.mod[which(sexes$Assignment=="MALE"), 3], na.rm=T)
  }

  #Run mosaic check per sample
  pMosaic <- as.data.frame(t(unlist(apply(sexes, 1, function(vals){
    if(any(is.na(vals))){
      pMos.X <- NA
      pMos.Y <- NA
    }else{
      #Get raw CNs
      CN.X <- dat.mod[which(dat.mod.rounded$sample_id==vals[1]), 2]
      CN.Y <- dat.mod[which(dat.mod.rounded$sample_id==vals[1]), 3]

      #Calculate p-value for chrX mosaicism
      if(CN.X>as.numeric(vals[2])){
        pMos.X <- pnorm(CN.X, mean=as.numeric(vals[2]), 
                        sd=sd.X*as.numeric(vals[2]), lower.tail=F)
      }else{
        pMos.X <- pnorm(CN.X, mean=as.numeric(vals[2]), 
                        sd=sd.X*as.numeric(vals[2]), lower.tail=T)
      }

      #Calculate p-value for chrY mosaicism
      if(CN.Y>as.numeric(vals[3])){
        pMos.Y <- pnorm(CN.Y, mean=as.numeric(vals[3]), 
                        sd=sd.Y*max(as.numeric(vals[3]), 1), lower.tail=F)
      }else{
        pMos.Y <- pnorm(CN.Y, mean=as.numeric(vals[3]), 
                        sd=sd.Y*max(as.numeric(vals[3]), 1), lower.tail=T)
      }
    }

    #Return p-values
    return(c(pMos.X, pMos.Y))
  }))))
  colnames(pMosaic) <- c("pMos.X", "pMos.Y")
  pMosaic$qMos.X <- p.adjust(pMosaic$pMos.X, method="fdr")
  pMosaic$qMos.Y <- p.adjust(pMosaic$pMos.Y, method="fdr")

  #Add mosaic check to sexes output
  sexes <- cbind(sexes, pMosaic)

  #####Plot sex assignments
  if(plot==T){
    #Prepare plot area
    par(mar=c(3.5, 3.5, 0.5, 0.5))
    plot(x=c(0, axLim), y=c(0, axLim), type="n", 
         xlab="", ylab="", xaxt="n", yaxt="n")

    #Add grid lines
    abline(h=0:axLim, v=0:axLim, lty=3, col="gray50")

    #Add sex karyotype labels behind each centroid
    sexCNs.df <- data.frame("X"=rep(0:axLim, axLim+1), 
                            "Y"=as.vector(unlist(sapply(0:axLim, rep, times=axLim+1))))
    sexCNs.df$karyo <- apply(sexCNs.df, 1, function(vals){
      return(paste(paste(rep("X", vals[1]), collapse=""), 
                   paste(rep("Y", vals[2]), collapse=""), 
                   sep=""))
    })
    text(x=sexCNs.df[, 1], y=sexCNs.df[, 2], 
         labels=sexCNs.df$karyo, 
         font=2, col="gray95", cex=0.8)

    #Assign colors for sex plotting
    if(mosaicThresh=="FDR"){
      if(length(sample.exclude)>0){
        colVect <- apply(sexes[-sample.exclude, c(4, 7:8)], 1, function(vals){
          if(vals[1] %in% sexAssign.df$label){
            return(sexAssign.df[which(sexAssign.df$label==vals[1]), ]$color)
          }else{
            return("#8F1336")
          }
        })
      }else{
        colVect <- apply(sexes[, c(4, 7:8)], 1, function(vals){
          if(vals[1] %in% sexAssign.df$label){
            return(sexAssign.df[which(sexAssign.df$label==vals[1]), ]$color)
          }else{
            return("#8F1336")
          }
        })
      }
    }else{
      if(length(sample.exclude)>0){
        colVect <- apply(sexes[-sample.exclude, c(4:6)], 1, function(vals){
          if(vals[1] %in% sexAssign.df$label){
            return(sexAssign.df[which(sexAssign.df$label==vals[1]), ]$color)
          }else{
            return("#8F1336")
          }
        })
      }else{
        colVect <- apply(sexes[, c(4:6)], 1, function(vals){
          if(vals[1] %in% sexAssign.df$label){
            return(sexAssign.df[which(sexAssign.df$label==vals[1]), ]$color)
          }else{
            return("#8F1336")
          }
        })
      }
    }

    #Plot points, highlighting specified sample
    isHighlightSample = dat.mod$sample_id == highlightSample
    pch_vector <- ifelse(isHighlightSample, 17, 19)
    cex_vector <- ifelse(isHighlightSample, 1.0, 0.5)
    points(dat.mod[, -1], pch=pch_vector, col=colVect, cex=cex_vector)
    abline(h = dat.mod[isHighlightSample, 3], lty=1, lwd=2, col="gray")
    abline(v = dat.mod[isHighlightSample, 2], lty=1, lwd=2, col="gray")

    #Add x-axis
    axis(1, at=0:axLim)
    mtext(1, line=2.2, text="chrX Copy Number")

    #Add y-axis
    axis(2, at=0:axLim, las=2)
    mtext(2, line=2.2, text="chrY Copy Number")

    #Add legend
    legendLabs <- apply(sexAssign.df[, 1:3], 1, function(vals){
      paste(vals[3], " (", 
            paste(rep("X", times=vals[1]), collapse=""), 
            paste(rep("Y", times=vals[2]), collapse=""), 
            ")", sep="")
    })
    legendLabs = c(legendLabs, "OTHER")
    colors = c(sexAssign.df$color, "#8F1336")
    legend("topright", bg="white",
           legend=legendLabs,
           pch=19, col=colors,
           cex=0.75, pt.cex=1)
    if (highlightSample != "" && highlightSample %in% dat.mod$sample_id) {
      legendLabs = c(legendLabs, highlightSample)
      colors = c(colors, "gray")
      legend("topleft", bg="white",
           legend=highlightSample,
           col="gray", lty=1, lwd=2, cex=1)
    }
  }

  #Return sex assignments
  return(sexes)
}

###############################################################
#####Helper function to plot distribution of samples per contig
###############################################################
boxplotsPerContig <- function(dat, exclude, genome.ploidy=2, contig.ploidy, 
                              contigLabels=paste("chr", c(1:22, "X", "Y"), sep=""), 
                              xmain="Chromosome", ymax=NULL, 
                              colorSignif=T, connect=F, boxes=T,
                              highlightSample=""){
  #Load library
  require(beeswarm)

  #Create plot color dataframe
  if(colorSignif==T){
    #Calculate p-values
    pvals <- testCNs(dat, FDR=T)

    #Iterate over pvalue matrix & compute plot color
    colMat <- sapply(2:ncol(pvals), function(col){
      sapply(1:nrow(pvals), function(row){
        if(is.na(pvals[row, col])){
          return(NA)
        }else{
          if(pvals[row, col]<0.05){
            if(dat[row, col]>contig.ploidy[col-1]){
              return("blue")
            }else{
              return("red")
            }
          }else{
            return("#838393")
          }
        }
      })
    })
  }else{
    colMat <- dat[, -1]
    colMat <- "#838393"
  }

  #Get max y-value
  if(is.null(ymax)){
    ymax <- max(4, max(dat[, -1], na.rm=T))
  }

  #Prepare plot area
  par(mar=c(3.5, 3.5, 0.5, 0.5))
  plot(x=c(0, ncol(dat)-1), y=c(0, ymax), 
       type="n", yaxt="n", xaxt="n", xaxs="i", ylab="", xlab="")

  #Add shading
  rect(xleft=seq(0, ncol(dat)+1, 2), xright=seq(1, ncol(dat)+2, 2), 
       ybottom=par("usr")[3], ytop=par("usr")[4], 
       border=NA, col="gray95")

  #Add gridlines
  abline(h=seq(0, ymax, 0.5), lty=2, col="gray85")
  abline(h=0:ymax, col="gray80")
  abline(h=c(0, genome.ploidy))

  #Iterate over contigs and plot all data
  if(connect==T){
    sapply(1:nrow(dat), function(i){
      vals <- dat[i, -1]
      cols <- colMat[i, ]
      sapply(1:length(vals), function(j){
        if(cols[j]!="#838393" & !is.na(cols[j])){
          if(j==1){
            segments(x0=j-0.5, x1=j+0.5, 
                     y0=as.numeric(vals[j]), 
                     y1=as.numeric(vals[j+1]), 
                     col="gray80", lwd=0.4)
          }else if(j==length(vals)){
            segments(x0=j-1.5, x1=j-0.5, 
                     y0=as.numeric(vals[j-1]), 
                     y1=as.numeric(vals[j]), 
                     col="gray80", lwd=0.4)
          }else{
            segments(x0=j-1.5, x1=j-0.5, 
                     y0=as.numeric(vals[j-1]), 
                     y1=as.numeric(vals[j]), 
                     col="gray80", lwd=0.4)
            segments(x0=j-0.5, x1=j+0.5, 
                     y0=as.numeric(vals[j]), 
                     y1=as.numeric(vals[j+1]), 
                     col="gray80", lwd=0.4)
          }
        }
      })
    })
  }
  sapply(1:(ncol(dat)-1), function(i){
    if(connect==F){
      #Jitter
      points(x=jitter(rep(i-0.5, times=nrow(dat)), amount=0.3), y=dat[, i+1], 
             pch=19, col=colMat[, i], cex=0.25)
      # #Swarm
      # beeswarm(dat[, i+1], add=T, at=i-0.5, method="swarm", 
      #          corral="wrap", corralWidth=0.6, 
      #          pch=19, pwcol=colMat[, i], cex=0.2)
    }else{
      points(x=rep(i-0.5, times=nrow(dat)), y=dat[, i+1], 
             pch=19, col=colMat[, i], cex=0.25)
    }
  })

  #Add boxplots
  if(boxes==T){
    if(!is.na(exclude)){
      boxplot(dat[, -c(1, exclude)], at=(1:(ncol(dat)-1))[-c(exclude-1)]-0.5, 
              add=T, outline=F, col=NA, lwd=0.75, lty=1, staplewex=0, 
              yaxt="n", xaxt="n", ylab="", xlab="")
    }else{
      boxplot(dat[, -1], at=(1:(ncol(dat)-1))-0.5, 
              add=T, outline=F, col=NA, lwd=0.75, lty=1, staplewex=0, 
              yaxt="n", xaxt="n", ylab="", xlab="")
    }
  }

  #Add single sample
  if (highlightSample != "" && highlightSample %in% dat[, 1]) {
    isHighlightSample = dat[, 1] == highlightSample
    lines((1:(ncol(dat)-1))-0.5, dat[isHighlightSample, -1], col="magenta")
    legend("topright", bg="white",
      legend=c(highlightSample),
      col=c("magenta"), lty=1)
  }

  #Add x-axis labels
  axis(1, at=(1:length(contigLabels))-0.5, tick=F, line=-0.8, las=2, labels=contigLabels)
  mtext(1, text=xmain, line=2.2)

  #Add y-axis labels
  axis(2, at=0:ymax, las=2)
  mtext(2, text="Estimated Copy Number", line=2.2)

  #Add text in top-left indicating number of samples
  nSamp <- nrow(dat)
  nSamp.hasNA <- length(unique(unlist(sapply(1:nrow(dat[, -1]), function(i){
    if(any(is.na(dat[i, -1]))){
      return(i)
    }else{
      return(NA)
    }
  }))))
  text(x=par("usr")[1], y=0.975*ymax, 
       labels=paste("N=", prettyNum(nSamp, big.mark=", "), " samples (", 
                    prettyNum(nSamp.hasNA-1, big.mark=", "), " incomplete)", 
                    sep=""), pos=4)

  #Add border cleanup
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4], 
       col=NA, border="black", lwd=2)
}

##########################
#####Rscript functionality
##########################
require(optparse)

#List of command-line options
option_list <- list(
  make_option(c("-O", "--OUTDIR"), type="character", default=NULL, 
              help="output directory [default: pwd]", 
              metavar="character"), 
  make_option(c("--noplot"), action="store_true", default=FALSE, 
              help="disable all visualization [default: %default]"), 
  make_option(c("-z", "--gzip"), action="store_true", default=FALSE, 
              help="gzip output files [default: %default]"), 
  make_option(c("-k", "--kmeans"), action="store_true", default=FALSE, 
              help="perform PCA & k-means based sample batching [default: %default]"), 
  make_option(c("-D", "--dimensions"), type="integer", default=8, 
              help="number of principal components to use in sample batching (requires -k) [default: %default]", 
              metavar="integer"), 
  make_option(c("--batchSize"), type="integer", default=100, 
              help="target number of samples per batch (requires -k) [default: %default]", 
              metavar="integer"), 
  make_option(c("--minBatch"), type="integer", default=60, 
              help="minimum number of samples per batch (requires -k) [default: %default]", 
              metavar="integer"), 
  make_option(c("--maxBatch"), type="integer", default=150, 
              help="maximum number of samples per batch (requires -k) [default: %default]", 
              metavar="integer"),
  make_option(c("--highlightSample"), default="",
              help="ID of sample to highlight in plots")
)

#Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog [options] coverage_matrix", 
                                option_list=option_list), 
                   positional_arguments=TRUE)

#Sanity-check arguments
if(length(args$args) != 1){
  stop("Must supply an input coverage matrix\n")
}

#Clean arguments & options
INFILE <- args$args[1]
OUTDIR <- args$options$OUTDIR
if(is.null(OUTDIR)){
  OUTDIR <- "./"
}
plot <- !(args$options$noplot)
gzip <- args$options$gzip
kmeans <- args$options$kmeans
nPCs <- args$options$dimensions
batch.ideal <- args$options$batchSize
batch.min <- args$options$minBatch
batch.max <- args$options$batchMax
highlightSample = args$options$highlightSample

# ##Jan 2020 dev parameters (on local machine)
# # INFILE <- "/Users/rlc/scratch/1KGP_2504_sub_batch_10_ploidy_matrix.bed.gz"
# INFILE <- "~/scratch/SFARI_Phase1_10_ploidy_matrix.bed.gz"
# plot <- T
# OUTDIR <- "~/scratch/WGDmodel_testing/"
# gzip <- T
# kmeans <- F
# nPCs <- 8
# batch.min <- 60
# batch.max <- 150
# batch.ideal <- 100

#Create OUTDIR if it doesn't already exist
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

#####PART 1: DATA PROCESSING#####
#Read, normalize, and clean coverage data
dat <- readMatrix(INFILE)
dat <- normalizeContigsPerSample(dat, exclude=c("X", "Y"), ploidy=2)
dat <- filterZeroBins(dat)
chr.dat <- medianPerContigPerSample(dat)
# chr.dat.norm <- normalizeContigsPerMatrix(chr.dat, scale.exclude=c("X", "Y"))

#####PART 2: DOSAGE-BASED BATCHING#####
#Only run if kmeans is optioned
if(kmeans==T){
  #Perform PCA on full matrix
  PCs <- binCovPCA(dat, exclude=c("X", "Y"), topPCs=nPCs)

  #Cluster samples based on dosage PCA
  #Note: tries this with seeds 1-100 (iterated sequentially) until first success
  #From what I can tell, this is necessary due to centroid initialization of k-means
  #This can probably be patched at a later date
  clust.PCs <- NULL
  seed <- 0
  while(is.null(clust.PCs) && seed<=100){
    seed <- seed+1
    set.seed(seed)
    try(
      if(plot==T){
        pdf(paste(OUTDIR, "/dosageClustering.scatter.pdf", sep=""), height=8, width=8*(7/6))
        clust.PCs <- clusterDosageProfiles(PCs$top[, 1:(nPCs+1)], plot=T)
        dev.off()
      }else{
        clust.PCs <- clusterDosageProfiles(PCs$top[, 1:(nPCs+1)], plot=F)
      }
      , silent=T)
    try(dev.off(), silent=T)
  }

  #Write out list of cluster assignments per sample
  clusterAssignments.out <- data.frame("sample_id"=names(clust.PCs$batch), 
                                       "batch"=as.vector(as.integer(clust.PCs$batch)))
  colnames(clusterAssignments.out)[1] <- "sample_id"
  write.table(clusterAssignments.out, 
              paste(OUTDIR, "/sample_batch_assignments.txt", sep=""), 
              col.names=T, row.names=F, sep="\t", quote=F)
  if(gzip==T){
    system(paste("gzip -f ", OUTDIR, "/sample_batch_assignments.txt", sep=""))
  }

  #Write out list of PCs per sample
  samplePCs.out <- PCs$top
  colnames(samplePCs.out)[1] <- "sample_id"
  write.table(samplePCs.out, paste(OUTDIR, "/sample_PCs.txt", sep=""), 
              col.names=T, row.names=F, sep="\t", quote=F)
  if(gzip==T){
    system(paste("gzip -f ", OUTDIR, "/sample_PCs.txt", sep=""))
  }

  #Write out list of PC coordinates per cluster center
  clusterCenters.out <- data.frame("batch"=1:nrow(clust.PCs$cluster.centers), 
                                   clust.PCs$cluster.centers)
  colnames(clusterCenters.out)[1] <- "#batch"
  write.table(clusterCenters.out, 
              paste(OUTDIR, "/batch_center_coordinates.txt", sep=""), 
              col.names=T, row.names=F, sep="\t", quote=F)
  if(gzip==T){
    system(paste("gzip -f ", OUTDIR, "/batch_center_coordinates.txt", sep=""))
  }
}

#####PART 3: SEX ASSIGNMENT#####
#Set sex assignment table
sexAssign.df <- data.frame("CN.X"=c(1, 2, 1, 3, 2, 1), 
                           "CN.Y"=c(1, 0, 0, 0, 1, 2), 
                           "label"=c("MALE", "FEMALE", "TURNER", 
                                     "TRIPLE X", "KLINEFELTER", "JACOBS"), 
                           "color"=c("#00BFF4", "#fd8eff", "#e02006", 
                                     "#7B2AB3", "#FF6A09", "#29840f"))

#Assign sexes
png(paste(OUTDIR, "/sex_assignments.png", sep=""), 
    height=1500, width=1500, res=300)
sexes <- assignSex(chr.dat, plot=T, sexAssign.df=sexAssign.df, highlightSample=highlightSample)
dev.off()

#Split data by chrX copy number (≥2 or <2) & evenly distribute NAs among M & F
sex.males <- as.vector(which(round(as.numeric(sexes$chrX.CN, 0))<2 & !is.na(sexes$chrX.CN)))
sex.females <- as.vector(which(round(as.numeric(sexes$chrX.CN, 0))>=2 & !is.na(sexes$chrX.CN)))
sex.NAs <- as.vector(which(is.na(sexes$chrX.CN)))
if(length(sex.NAs)==0){
  chr.dat.males <- chr.dat[c(sex.males), ]
  chr.dat.females <- chr.dat[c(sex.females), ]
}else{
  if(length(sex.NAs)>1){
    sex.NAs.first <- sex.NAs[1:floor((length(sex.NAs)/2))]
    sex.NAs.second <- sex.NAs[(floor((length(sex.NAs)/2))+1):length(sex.NAs)]
  }else if(length(sex.NAs)==1){
    sex.NAs.first <- sex.NAs
    sex.NAs.second <- NULL
  }
  chr.dat.males <- chr.dat[c(sex.males, sex.NAs.first), ]
  chr.dat.females <- chr.dat[c(sex.females, sex.NAs.second), ]
}

#####PART 4: CNA/CNV SCREEN#####
#Plot CN per contig - Males
png(paste(OUTDIR, "/estimated_CN_per_contig.chrX_lessThan_2copies.with_contours.png", sep=""), 
    height=1250, width=2500, res=300)
boxplotsPerContig(chr.dat.males, exclude=NA, contig.ploidy=c(rep(2, 22), 1, 1), connect=T, highlightSample=highlightSample)
dev.off()
png(paste(OUTDIR, "/estimated_CN_per_contig.chrX_lessThan_2copies.no_contours.png", sep=""), 
    height=1250, width=2500, res=300)
boxplotsPerContig(chr.dat.males, exclude=NA, contig.ploidy=c(rep(2, 22), 1, 1), connect=F, highlightSample=highlightSample)
dev.off()

#Plot CN per contig - Females
png(paste(OUTDIR, "/estimated_CN_per_contig.chrX_atLeast_2copies.with_contours.png", sep=""), 
    height=1250, width=2500, res=300)
boxplotsPerContig(chr.dat.females, exclude=NA, contig.ploidy=c(rep(2, 22), 2, 0), connect=T, highlightSample=highlightSample)
dev.off()
png(paste(OUTDIR, "/estimated_CN_per_contig.chrX_atLeast_2copies.no_contours.png", sep=""), 
    height=1250, width=2500, res=300)
boxplotsPerContig(chr.dat.females, exclude=NA, contig.ploidy=c(rep(2, 22), 2, 0), connect=F, highlightSample=highlightSample)
dev.off()

#Plot CN per contig - all samples
png(paste(OUTDIR, "/estimated_CN_per_contig.all_samples.with_contours.png", sep=""), 
    height=1250, width=2500, res=300)
boxplotsPerContig(chr.dat, exclude=NA, contig.ploidy=c(rep(2, 22), 2, 0), connect=T, boxes=F, highlightSample=highlightSample)
dev.off()
png(paste(OUTDIR, "/estimated_CN_per_contig.all_samples.no_contours.png", sep=""), 
    height=1250, width=2500, res=300)
boxplotsPerContig(chr.dat, exclude=NA, contig.ploidy=c(rep(2, 22), 2, 0), connect=F, boxes=F, highlightSample=highlightSample)
dev.off()

#Write table of sexes
sexes <- sexes[match(colnames(dat)[-c(1:3)], sexes$sample_id), ]
colnames(sexes)[1] <- "sample_id"
write.table(sexes, paste(OUTDIR, "/sample_sex_assignments.txt", sep=""), 
            col.names=T, row.names=F, sep="\t", quote=F)
if(gzip==T){
  system(paste("gzip -f ", OUTDIR, "/sample_sex_assignments.txt", sep=""))
}

#Generate p-values, q-values, and rounded CNs for males/females for whole chromosomes
males.p <- testCNs(chr.dat.males, FDR=F)
males.q <- testCNs(chr.dat.males, FDR=T)
males.CN <- chr.dat.males
males.CN[, -1] <- apply(males.CN[, -1], 2, round, digits=2)
females.p <- testCNs(chr.dat.females, FDR=F)
females.q <- testCNs(chr.dat.females, FDR=T)
females.CN <- chr.dat.females
females.CN[, -1] <- apply(females.CN[, -1], 2, round, digits=2)

#Merge male/female p-values and rounded CNs
merged.p <- rbind(males.p, females.p)
merged.p <- merged.p[match(colnames(dat)[-c(1:3)], merged.p$sample_id), ]
colnames(merged.p) <- c("sample_id", paste("chr", c(1:22, "X", "Y"), "_pValue", sep=""))
merged.q <- rbind(males.q, females.q)
merged.q <- merged.q[match(colnames(dat)[-c(1:3)], merged.q$sample_id), ]
colnames(merged.q) <- c("sample_id", paste("chr", c(1:22, "X", "Y"), "_qValue", sep=""))
merged.CN <- rbind(males.CN, females.CN)
merged.CN <- merged.CN[match(colnames(dat)[-c(1:3)], merged.CN$sample_id), ]
colnames(merged.CN) <- c("sample_id", paste("chr", c(1:22, "X", "Y"), "_CopyNumber", sep=""))

#Write merged p-values
write.table(merged.p, paste(OUTDIR, "/CNA_pValues.txt", sep=""), 
            col.names=T, row.names=F, sep="\t", quote=F)
if(gzip==T){
  system(paste("gzip -f ", OUTDIR, "/CNA_pValues.txt", sep=""))
}

#Write merged q-values
write.table(merged.q, paste(OUTDIR, "/CNA_qValues.txt", sep=""), 
            col.names=T, row.names=F, sep="\t", quote=F)
if(gzip==T){
  system(paste("gzip -f ", OUTDIR, "/CNA_qValues.txt", sep=""))
}

#Write merged copy number estimates
write.table(merged.CN, paste(OUTDIR, "/estimated_copy_numbers.txt", sep=""), 
            col.names=T, row.names=F, sep="\t", quote=F)
if(gzip==T){
  system(paste("gzip -f ", OUTDIR, "/estimated_copy_numbers.txt", sep=""))
}

#Plot binwise CN estimates per autosome -- all samples
binwise.dat <- data.frame(colnames(dat)[-c(1:3)], 
                          t(dat[, -c(1:3)]))
bin.IDs <- as.character(apply(dat[1:3], 1, paste, collapse="_", sep=""))
bin.IDs <- as.character(sapply(bin.IDs, function(ID){gsub(" ", "", ID, fixed=T)}))
colnames(binwise.dat) <- c("sample_id", bin.IDs)
sapply(setdiff(unique(dat[, 1]), c("X", "Y")), function(contig){
  png(paste(OUTDIR, "/estimated_CN_per_bin.all_samples.", contig, ".png", sep=""), 
      height=1250, width=2500, res=300)
  plot.dat <- binwise.dat[, c(1, which(dat[, 1]==contig)+1)]
  boxplotsPerContig(plot.dat, exclude=NA, contig.ploidy=rep(2, ncol(plot.dat)-1), connect=T, 
                    xmain=paste(contig, "Position (Binned)"), ymax=5, contigLabels=NA, highlightSample=highlightSample)
  dev.off()
})

#Plot binwise CN estimates per sex chromosome -- all samples
binwise.dat.males <- data.frame(chr.dat.males$sample_id, 
                                t(dat[, which(colnames(dat) %in% chr.dat.males$sample_id)]))
colnames(binwise.dat.males) <- c("sample_id", bin.IDs)
binwise.dat.females <- data.frame(chr.dat.females$sample_id, 
                                t(dat[, which(colnames(dat) %in% chr.dat.females$sample_id)]))
colnames(binwise.dat.females) <- c("sample_id", bin.IDs)
sapply(intersect(unique(dat[, 1]), c("X", "Y")), function(contig){
  #Males
  png(paste(OUTDIR, "/estimated_CN_per_bin.chrX_lessThan_2copies.", contig, ".png", sep=""), 
      height=1250, width=2500, res=300)
  plot.dat <- binwise.dat.males[, c(1, which(dat[, 1]==contig)+1)]
  boxplotsPerContig(plot.dat, exclude=NA, contig.ploidy=rep(2, ncol(plot.dat)-1), connect=T, 
                    xmain=paste(contig, "Position (Binned)"), ymax=5, contigLabels=NA, highlightSample=highlightSample)
  dev.off()
  #Females
  png(paste(OUTDIR, "/estimated_CN_per_bin.chrX_atLeast_2copies.", contig, ".png", sep=""), 
      height=1250, width=2500, res=300)
  plot.dat <- binwise.dat.females[, c(1, which(dat[, 1]==contig)+1)]
  boxplotsPerContig(plot.dat, exclude=NA, contig.ploidy=rep(2, ncol(plot.dat)-1), connect=T, 
                    xmain=paste(contig, "Position (Binned)"), ymax=5, contigLabels=NA, highlightSample=highlightSample)
  dev.off()
})

#Generate p-values, q-values, and rounded CNs for males/females for per-bin
males.binwise.p <- testCNs(binwise.dat.males, FDR=F)
males.binwise.q <- testCNs(binwise.dat.males, FDR=T)
males.binwise.CN <- binwise.dat.males
males.binwise.CN[, -1] <- apply(males.binwise.CN[, -1], 2, round, digits=2)
females.binwise.p <- testCNs(binwise.dat.females, FDR=F)
females.binwise.q <- testCNs(binwise.dat.females, FDR=T)
females.binwise.CN <- binwise.dat.females
females.binwise.CN[, -1] <- apply(females.binwise.CN[, -1], 2, round, digits=2)

#Merge male/female p-values and rounded CNs
merged.p <- rbind(males.binwise.p, females.binwise.p)
merged.p <- merged.p[match(colnames(dat)[-c(1:3)], merged.p$sample_id), ]
merged.p <- t(merged.p)
merged.p <- cbind(dat[, 1:3], merged.p[-1, ])
colnames(merged.p)[1] <- c("#Chr")
merged.q <- rbind(males.binwise.q, females.binwise.q)
merged.q <- merged.q[match(colnames(dat)[-c(1:3)], merged.q$sample_id), ]
merged.q <- t(merged.q)
merged.q <- cbind(dat[, 1:3], merged.q[-1, ])
colnames(merged.q)[1] <- c("#Chr")
merged.CN <- rbind(males.binwise.CN, females.binwise.CN)
merged.CN <- merged.CN[match(colnames(dat)[-c(1:3)], merged.CN$sample_id), ]
merged.CN <- t(merged.CN)
merged.CN <- cbind(dat[, 1:3], merged.CN[-1, ])
colnames(merged.CN)[1] <- c("#Chr")

#Write merged p-values
write.table(merged.p, paste(OUTDIR, "/binwise_CNV_pValues.bed", sep=""), 
            col.names=T, row.names=F, sep="\t", quote=F)
if(gzip==T){
  system(paste("gzip -f ", OUTDIR, "/binwise_CNV_pValues.bed", sep=""))
}

#Write merged q-values
write.table(merged.q, paste(OUTDIR, "/binwise_CNV_qValues.bed", sep=""), 
            col.names=T, row.names=F, sep="\t", quote=F)
if(gzip==T){
  system(paste("gzip -f ", OUTDIR, "/binwise_CNV_qValues.bed", sep=""))
}

#Write merged copy number estimates
write.table(merged.CN, paste(OUTDIR, "/binwise_estimated_copy_numbers.bed", sep=""), 
            col.names=T, row.names=F, sep="\t", quote=F)
if(gzip==T){
  system(paste("gzip -f ", OUTDIR, "/binwise_estimated_copy_numbers.bed", sep=""))
}

