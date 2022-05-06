#' Manually combine multiple runs of SpaRTaCo
#'
#' This function combines multiple runs of SpaRTaCo with the same values of K and R.

#' @import ggplot2
#' @export
#'
#' @param x either a list of `spartaco` objects or a vector containing the names of the files to be combined;
#' @param search.dir the directory path where the file names given by `x` are searched. If `NULL`, the file paths must be specified as part of the file names is `x`;
#' @param compute.uncertainty if `TRUE` (default), it computes the clustering uncertainty of the rows and of the columns.
#'
#' @return an object of class `spartaco` with the parameter estimates and the clustering labels given by the best fit across the ones in `x`. If `compute.uncertainty == T`, the clustering uncertainty measures are returned.
#'
#' @details If each element of `x` is an object of class `spartaco`, then the function combines them into a unique spartaco object.
#' If the `spartaco` objects are stored into different files, the vector `x` receives the names of the files,
#' and `search.dir` is the path of the directory where the files are stored.

#' @references
#' Sottosanti, A. and Risso, D. (2021+) Co-clustering of Spatially Resolved Transcriptomic Data [(preprint)](https://arxiv.org/abs/2110.04872)

#' @seealso [spartaco::spartaco()]
#'
#' @examples
#' library(spartaco)
#'
#' # The following code illustruates how to combine three output of the spartaco function saved in the files output1.Rdata, output2.Rdata and output3.Rdata.
#' # The directory containing the files is specified as follows:
#' directory <- "..."
#' x <- c("output1.Rdata", "output2.Rdata", "output3.Rdata")
#' CombineSpartaco(x = x, search.dir = directory)
#'
#' # In alternative, one could run:
#' x <- paste(directory, x, sep = "/")
#' CombineSpartaco(x = x)

CombineSpartaco <- function(x = NULL, #KR = NULL,
                           search.dir = NULL, compute.uncertainty = T){
    loadRData <- function(fileName){
        #loads an RData file, and returns it
        load(fileName)
        get(ls()[ls() != "fileName"])
    }

    if(is.list(x)) results <- x else {
        results <- list()
        for(i in 1:length(x)){
            file.name <- x[i]
            if(!is.null(search.dir)){
                if(substr(search.dir, nchar(search.dir),nchar(search.dir)) == "/")
                    file.name <- paste(search.dir, file.name, sep="") else
                        file.name <- paste(search.dir, file.name, sep="/")
            }
            results[[i]] <- loadRData(file.name)
        }
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
            final$stn.ratio <- final$tau/final$xi
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
    if(compute.uncertainty){
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

    # execution time
    final$exec.time <- numeric(length(x))
    sapply(1:length(x), function(i)
        if(is.null(results[[i]]$exec.time)) final$exec.time[i] <<- NA else
            final$exec.time[i] <<- results[[i]]$exec.time)

    ChainPlot <- data.frame(
        ranges = as.vector(unlist(likelihoods)),
        iterations = as.vector(unlist(sapply(1:length(x), function(j) 1:length(likelihoods[[j]])))),
        chains = as.factor(as.vector(unlist(sapply(1:length(x), function(j) rep(j, length(likelihoods[[j]]))))))
    )
    final$ll.plot <- ggplot(ChainPlot, aes(iterations, ranges, col = chains))+geom_point()+
        theme_classic()+theme(legend.position = "bottom")+
        labs(x = "iteration", y = "classification log-likelihood", color = "run")
    plot(final$ll.plot)
       # plot(1,1,
        #     xlim = c(1, max(sapply(1:length(likelihoods), function(l) length(likelihoods[[l]])))),
        #     ylim = c(min(ranges), max(ranges)),
        #     xlab = "Iteration", ylab = expression(logL))
        #for(i in 1:length(likelihoods))
        #    points(1:length(likelihoods[[i]]), likelihoods[[i]], col = i, pch = i)

    class(final) <- "spartaco"
    return(final)
}
