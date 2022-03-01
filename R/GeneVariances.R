#' Compute the gene-specific variances
#'
#' This function returns the summary statistics for the random variances estimated through SpaRTaCo within the `K*R` blocks.
#'
#' @import SpatialExperiment
#' @import invgamma
#' @import coda
#' @export
#'
#' @param results The output of a `spartaco` function;
#' @param gene.names A vector of length `nrow(results$x)`, containing the names of the genes. If `NULL`, the function takes `gene.names = row.names(results$x)`.
#' @param HPD.interval A list containing the options for computing the random parameters high posterior density intervals of level `HPD.interval$level` (default `0.95`).
#' The intervals are computed numerically via MCMC; the number of draws used to compute the intervals is given by `HPD.interval$level` (default is `10^4`).
#' If `NULL`, the intervals are not computed.
#'
#' @return A list of five elements:
#' - `Expectation`: a data frame of dimension `nrow(results$x)` x `R` containing the expected values of the random effects within the `R` spot clusters;
#' - `Variance`: a data frame of dimension `nrow(results$x)` x `R` containing the variance of the random effects within the `R` spot clusters;
#' - `HPD.left`: a data frame of dimension `nrow(results$x)` x `R` containing the left HPD interval of the random effects within the `R` spot clusters (only if `!is.null(HPD.interval)`);
#' - `HPD.right` a data frame of dimension `nrow(results$x)` x `R` containing the right HPD interval of the random effects within the `R` spot clusters (only if `!is.null(HPD.interval)`);
#' - `Cs`: the row clustering labels (they are the same contained into `results$Cs`).
#'
#' More details on how to access the random parameter estimates are given in the Examples.
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

GeneVariances <- function(results, gene.names = NULL, HPD.interval = list(level = .95, R = 10^4)){
    x <- results$x
    coordinates <- results$coord
    Mu <- results$mu
    Tau <- results$tau
    Xi <- results$xi
    Alpha <- results$alpha
    Beta <- results$beta
    Phi <- results$phi
    Cs <- results$Cs
    Ds <- results$Ds
    if(is.null(gene.names)) gene.names <- row.names(x)
    if(is.null(gene.names)) gene.names <- 1:nrow(x)
    expDist <- as.matrix(exp(-as.matrix(dist(coordinates))))
    K <- ifelse(is.vector(Mu), 1, nrow(Mu))
    R <- length(Phi)
    Alpha.post <- Beta.post <- Expected.post <- Variance.post <- data.frame(matrix(0, nrow(x), R))
    rownames(Alpha.post) <- rownames(Beta.post) <- rownames(Expected.post) <-
        rownames(Variance.post) <- gene.names
    colnames(Alpha.post) <- colnames(Beta.post) <- colnames(Expected.post) <-
        colnames(Variance.post) <- paste("r =",1:R)
    if(!is.null(HPD.interval)){
        Inter.low <- Inter.up <- data.frame(matrix(0, nrow(x), R))
        rownames(Inter.low) <- rownames(Inter.up) <- gene.names
        colnames(Inter.low) <- colnames(Inter.up) <- paste("r =",1:R)
    }
    for(k in sort(unique(Cs))){
        for(r in 1:R){
            eig.r <- eigen(expDist[Ds == r, Ds == r]^(1/Phi[r]))
            Block1 <- (x[Cs == k, Ds == r] - Mu[k, r]) %*% eig.r$vec
            invD <- (1/(Tau[k, r]*eig.r$val + Xi[k, r]))
            Alpha.post[Cs == k,r] <- sum(Ds == r)/2 + Alpha[k, r]
            Beta.post[Cs == k,r] <- (Block1 * Block1) %*% invD/2 + Beta[k, r]
            Expected.post[Cs == k,r] <- ifelse(Alpha.post[Cs == k,r]>1, Beta.post[Cs == k,r]/(Alpha.post[Cs == k,r]-1), Inf)
            Variance.post[Cs == k,r] <- ifelse(Alpha.post[Cs == k,r]>2, Beta.post[Cs == k,r]^2/((Alpha.post[Cs == k,r]-1)^2 * (Alpha.post[Cs == k,r]-2)), Inf)
            if(!is.null(HPD.interval)){
                for(i in 1:sum(Cs == k)){
                    gen <- rinvgamma(HPD.interval$R, shape = Alpha.post[Cs == k,r][i], rate = Beta.post[Cs == k,r][i])
                    interv <- HPDinterval(obj = mcmc(gen), prob = HPD.interval$level)
                    Inter.low[Cs == k,r][i] <- interv[1]
                    Inter.up[Cs == k,r][i] <- interv[2]
                }
            }
        }
    }
    if(is.null(HPD.interval)) Inter.low <- Inter.up <- NULL
    return(list(Expectation = Expected.post,
                Variance = Variance.post,
                HPD.left = Inter.low,
                HPD.right = Inter.up,
                Cs = Cs))
}
