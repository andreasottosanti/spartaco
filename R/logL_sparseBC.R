#' Log-likelihood and ICL of the Sparse Biclustering model
#'
#' This function computes the log-likelihood and the ICL of a sparse biclustering model (see *Tan and Witten, 2014*). The ICL is computed as described by *Sottosanti and Risso (2021+)*.
#'
#' @export
#'
#' @param x the data matrix;
#' @param lambda the penalization parameter;
#' @param results an object of class `sparseBC`;
#' @param Cs the row clustering labels (only if `results == NULL`);
#' @param Ds the column clustering labels (only if `results == NULL`);
#' @param mu the `K` x `R` matrix containing the co-cluster centroids (only if `results == NULL`);
#'
#' @return A list containing the log-likelihood, the ICL and the ML estimate of the scale parameter, that it is same for every block.
#'
#' @details If an object of class `sparseBC` is passed to `results`, the values of `Cs`, `Ds` and `mu`. The penalization `lambda` needs always to be specified, as it cannot be extracted from the `sparseBC` object.
#'
#' @seealso [sparseBC::sparseBC()]
#'
#' @references
#' Sottosanti, A. and Risso, D. (2021+) Co-clustering of Spatially Resolved Transcriptomic Data [(preprint)](https://arxiv.org/abs/2110.04872).
#'
#' Tan, K. M., and Witten, D. M. (2014) Sparse biclustering of transposable data. *Journal of Computational and Graphical Statistics* **23(4)** 985-1008.
#'
#' @examples
#' library(spartaco)
#' library(sparseBC)
#'
#' x <- matrix(runif(n*p), n, p)
#' results <- sparseBC(x, k = k, r = r, lambda = lambda)
#'
#' The following lines are equivalent:
#' logL.sparseBC(x = x, results = results, lambda = lambda)
#' logL.sparseBC(x = x, Cs = results$Cs, Ds = results$Ds, mu = results$Mus, lambda = lambda)
#'

logL.sparseBC <- function(x, lambda, results = NULL, Cs = NULL, Ds = NULL, mu = NULL){
    if(!is.null(results)){
        if(attr(results, "class") == "sparseBC"){
            Cs <- results$Cs
            Ds <- results$Ds
            mu <- results$Mu
        }
    }
    nk <- table(Cs)
    pr <- table(Ds)
    K <- length(nk)
    R <- length(pr)
    sigma2.num <- sigma2.den <- matrix(0, K, R)
    sapply(1:K, function(k){
        sapply(1:R, function(r){
            sigma2.num[k,r] <<- sum((x[Cs == k, Ds == r]-mu[k,r])^2)
            sigma2.den[k,r] <<- nk[k]*pr[r]
        })
    })
    sigma2 <- sum(sigma2.num)/sum(sigma2.den)
    logL <- matrix(0, K, R)
    sapply(1:K, function(k){
        sapply(1:R, function(r){
            logL[k,r] <<- dmvn(x[Cs == k, Ds == r], mu = mu[k,r], sigma = sigma2*diag(1, nrow = nk[k]), delta = diag(1, nrow = pr[r]))
        })
    })
    return(list(
        logL = sum(logL)-lambda*sum(abs(mu)),
        ICL = (sum(logL)-lambda*sum(abs(mu)))-nrow(x)*log(K)-ncol(x)*log(R)-(K*R+1)/2*log(prod(dim(x))),
        sigma2 = sigma2))
}
