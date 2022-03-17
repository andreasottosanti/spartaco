CombineResults <- function(results = NULL, KR = NULL, search.dir = NULL, display = F, compute.discrepancy = F){
    loadRData <- function(fileName){
        #loads an RData file, and returns it
        load(fileName)
        get(ls()[ls() != "fileName"])
    }
    if(is.null(results)){
        if(is.null(search.dir)) search.dir <- getwd()
        dat <- dir(search.dir)
        dat <- dat[grepl(paste("K",KR[1],"_R",KR[2], sep = ""), dat)]
        cat(paste("Found", length(dat),"files in the current directory\n"))
        results <- list()
        for(i in 1:length(dat)) results[[i]] <- loadRData(paste(search.dir,dat[i],sep="/"))
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
    if(compute.discrepancy){
        DS <- sapply(1:length(results), function(i) results[[i]]$Ds)
        final$CERs <- diag(0,length(results),length(results))
        sapply(1:(length(results)-1), function(i){
            sapply(i:length(results), function(j){
                final$CERs[i,j] <<- final$CERs[j,i] <<- CER(DS[,i], DS[,j])
            })
        })
    }
    final$max.logL <- maxi
    final$x <- results[[i]]$x
    final$coordinates <- results[[i]]$coordinates
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
