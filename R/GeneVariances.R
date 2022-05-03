#' Compute the gene-specific variances
#'
#' This function returns the summary statistics for the random variances estimated through SpaRTaCo within the `K*R` blocks.
#'
# #' @import SpatialExperiment
#' @import invgamma
#' @import coda
#' @export
#'
#' @param x an object of class `spartaco`;
#' @param gene.names a vector of length `nrow(x$x)`, containing the gene names. If `NULL`, the function takes `gene.names = row.names(x$x)`.
#' @param MCMC.interval a list of length 2 containing the options for computing the high posterior density intervals of the gene specific variances. The  level is given by `MCMC.interval$level` (default `0.95`).
#' Intervals are computed via MCMC; the number of simulated draws must is passed to `MCMC.interval$level` (default is `10^4`).
#' If `NULL`, the intervals are computed as mean + c(-1.96,1.96)*sd, and the function returns a warning.
#'
#' @return An object of class `spartaco.genes` containing five elements:
#' - `Expectation`: a data frame of dimension `nrow(x$x)` x `R` containing the expected values of the random effects within the `R` spot clusters;
#' - `Variance`: a data frame of dimension `nrow(x$x)` x `R` containing the variance of the random effects within the `R` spot clusters;
#' - `HPD.left`: a data frame of dimension `nrow(x$x)` x `R` containing the left HPD interval of the random effects within the `R` spot clusters;
#' - `HPD.right` a data frame of dimension `nrow(x$x)` x `R` containing the right HPD interval of the random effects within the `R` spot clusters;
#' - `Cs`: the row clustering labels (the same contained into the object `x`).
#'
#' More details on how to access the random parameter estimates are given in the **Examples**.
#'
#' @references
#'  Sottosanti, A. and Risso, D. (2021+) Co-clustering of Spatially Resolved Transcriptomic Data [(preprint)](https://arxiv.org/abs/2110.04872)
#'
#' @examples
#' library(spartaco)
#' # Not run:
#' # x <- matrix(runif(n*p), n, p)
#' # coordinates <- matrix(runif(2*p), p, 2)
#' # output <- spartaco(x = x, coordinates = coordinates, K = K, R = R)
#'
#' sigma2 <- GeneVariances(output)
#'
#' # To access the expectations of genes within the kr-th block:
#' sigma2$Expectation[sigma2$Cs == k,r]
#'
#' # To access the variances of the genes within the kr-th block:
#' sigma2$Variance[sigma2$Cs == k,r]
#'
#' # To access the HPD interval limits of the genes within the kr-th block:
#' sigma2$HPD.left[sigma2$Cs == k,r]
#' sigma2$HPD.right[sigma2$Cs == k,r]

GeneVariances <- function(x, gene.names = NULL, MCMC.interval = list(level = .95, R = 10^4)){
    if(class(x) != "spartaco") stop("x is not a spartaco object")
    X <- x$x
    coordinates <- x$coord
    Mu <- x$mu
    Tau <- x$tau
    Xi <- x$xi
    Alpha <- x$alpha
    Beta <- x$beta
    Phi <- x$phi
    Cs <- x$Cs
    Ds <- x$Ds
    if(is.null(gene.names)) gene.names <- row.names(X)
    if(is.null(gene.names)) gene.names <- 1:nrow(X)
    expDist <- as.matrix(exp(-as.matrix(dist(coordinates))))
    K <- ifelse(is.vector(Mu), 1, nrow(Mu))
    R <- length(Phi)
    Alpha.post <- Beta.post <- Expected.post <- Variance.post <- data.frame(matrix(0, nrow(X), R))
    rownames(Alpha.post) <- rownames(Beta.post) <- rownames(Expected.post) <-
        rownames(Variance.post) <- gene.names
    colnames(Alpha.post) <- colnames(Beta.post) <- colnames(Expected.post) <-
        colnames(Variance.post) <- paste("r =",1:R)

    Inter.low <- Inter.up <- data.frame(matrix(0, nrow(X), R))
    rownames(Inter.low) <- rownames(Inter.up) <- gene.names
    colnames(Inter.low) <- colnames(Inter.up) <- paste("r =",1:R)

    for(k in sort(unique(Cs))){
        for(r in 1:R){
            eig.r <- eigen(expDist[Ds == r, Ds == r]^(1/Phi[r]))
            Block1 <- (X[Cs == k, Ds == r] - Mu[k, r]) %*% eig.r$vec
            invD <- (1/(Tau[k, r]*eig.r$val + Xi[k, r]))
            Alpha.post[Cs == k,r] <- sum(Ds == r)/2 + Alpha[k, r]
            Beta.post[Cs == k,r] <- (Block1 * Block1) %*% invD/2 + Beta[k, r]
            Expected.post[Cs == k,r] <- ifelse(Alpha.post[Cs == k,r]>1, Beta.post[Cs == k,r]/(Alpha.post[Cs == k,r]-1), Inf)
            Variance.post[Cs == k,r] <- ifelse(Alpha.post[Cs == k,r]>2, Beta.post[Cs == k,r]^2/((Alpha.post[Cs == k,r]-1)^2 * (Alpha.post[Cs == k,r]-2)), Inf)
            if(!is.null(MCMC.interval)){
                for(i in 1:sum(Cs == k)){
                    gen <- rinvgamma(MCMC.interval$R, shape = Alpha.post[Cs == k,r][i], rate = Beta.post[Cs == k,r][i])
                    interv <- HPDinterval(obj = mcmc(gen), prob = MCMC.interval$level)
                    Inter.low[Cs == k,r][i] <- interv[1]
                    Inter.up[Cs == k,r][i] <- interv[2]
                }
            } else {
                Inter.low[,r] <- Expected.post[,r] - 2*sqrt(Variance.post[,r])
                Inter.up[,r] <- Expected.post[,r] + 2*sqrt(Variance.post[,r])
            }
        }
    }
    if(is.null(MCMC.interval)) warning("the HPD intervals returned are approximated as mean + c(-1.96,1.96)*sd")
    output <- list(Expectation = Expected.post,
                   Variance = Variance.post,
                   HPD.left = Inter.low,
                   HPD.right = Inter.up,
                   Cs = Cs)
    class(output) <- "spartaco.genes"
    return(output)
}
