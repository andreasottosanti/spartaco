CombineResults <- function(x = NULL, KR = NULL, search.dir = NULL, display = F, compute.discrepancy = F){
    loadRData <- function(fileName){
        #loads an RData file, and returns it
        load(fileName)
        get(ls()[ls() != "fileName"])
    }
    if(is.null(x)){
        if(is.null(search.dir)) search.dir <- getwd()
        dat <- dir(search.dir)
        dat <- dat[grepl(paste("K",KR[1],"_R",KR[2], sep = ""), dat)]
        cat(paste("Found", length(dat),"files in the current directory\n"))
        results <- list()
        for(i in 1:length(dat)) results[[i]] <- loadRData(paste(search.dir,dat[i],sep="/"))
    } else {
        results <- list()
        for(i in 1:length(dat)) results[[i]] <- loadRData(x[i])
    }
    likelihoods <- list()
    maxi <- rep(-Inf, length(results))
    final <- list()
    for(i in 1:length(results)){
        likelihoods[[i]] <- results[[i]]$logL
        maxi[i] <- max(likelihoods[[i]])
        if(which.max(maxi) == i){
            final$mu <- results[[i]]$mu
            final$tau <- results[[i]]$tau
            final$xi <- results[[i]]$xi
            final$alpha <- results[[i]]$alpha
            final$beta <- results[[i]]$beta
            final$phi <- results[[i]]$phi
            final$Cs <- results[[i]]$Cs
            final$Ds <- results[[i]]$Ds
            final$ICL <- results[[i]]$ICL
        }
    }
    final$max.logL <- maxi
    final$x <- results[[i]]$x
    final$coordinates <- results[[i]]$coordinates
    #final$CS <- sapply(1:length(results), function(i) results[[i]]$Cs)
    #final$DS <- sapply(1:length(results), function(i) results[[i]]$Ds)
    if(compute.discrepancy){
        best.j <- which.max(final$max.logL)
        j.to.invest <- setdiff(1:length(final$max.logL), best.j)
        final$cluster.discr <- list()
        # evaluate the discrepancy across the rows
        cat("Computing the row discrepancy...\n")
        CERs <- matrix(0, nrow(final$mu), length(j.to.invest))
        sapply(1:nrow(final$mu), function(k){
            reference <- as.numeric(final$Cs == k)
            sapply(1:length(j.to.invest), function(j){
                classif <- table(final$Cs, results[[j.to.invest[j]]]$Cs)[k,]
                k.j <- as.numeric(attr(classif, "names"))[which.max(classif)]
                comparison <- as.numeric(results[[j.to.invest[j]]]$Cs == k.j)
                CERs[k,j] <<- CER(reference = reference, estimate = comparison)
            })
        })
        w <- 1/(final$max.logL[best.j]-final$max.logL[-best.j])
        final$cluster.discr$rows <- as.vector((CERs %*% w)/sum(w))

        # evaluate the discrepancy across the columns
        cat("Computing the column discrepancy...\n")
        CERs <- matrix(0, ncol(final$mu), length(j.to.invest))
        sapply(1:ncol(final$mu), function(r){
            reference <- as.numeric(final$Ds == r)
            sapply(1:length(j.to.invest), function(j){
                classif <- table(final$Ds, results[[j.to.invest[j]]]$Ds)[r,]
                r.j <- as.numeric(attr(classif, "names"))[which.max(classif)]
                comparison <- as.numeric(results[[j.to.invest[j]]]$Ds == r.j)
                CERs[r,j] <<- CER(reference = reference, estimate = comparison)
            })
        })
        w <- 1/(final$max.logL[best.j]-final$max.logL[-best.j])
        final$cluster.discr$columns <- as.vector((CERs %*% w)/sum(w))
    }

    if(display){
        ranges <- as.vector(unlist(likelihoods))
        plot(1,1,
             xlim = c(1, max(sapply(1:length(likelihoods), function(l) length(likelihoods[[l]])))),
             ylim = c(min(ranges), max(ranges)),
             xlab = "Iteration", ylab = expression(logL))
        for(i in 1:length(likelihoods))
            points(1:length(likelihoods[[i]]), likelihoods[[i]], col = i, pch = i)
    }
    class(final) <- "spartaco"
    return(final)
}
