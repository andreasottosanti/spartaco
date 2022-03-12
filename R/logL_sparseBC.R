#' Log-likelihood and ICL of the Sparse Biclustering model
#'
#' This function computes the log-likelihood and the ICL of the Sparse Biclustering model of Tan & Witten (2014).
#'
#' @export
#'
#' @param x the data matrix;
#' @param results an object of class `sparseBC`;
#' @param Cs the row clustering labels (only if `results = NULL`);
#' @param Ds the column clustering labels (only if `results = NULL`);
#' @param mu the `K`x`R` matrix containing the co-cluster centroids (only if `results = NULL`);
#' @param lambda the penalization parameter (only if `results = NULL`);
#'
#' @return A list containing the log-likelihood, the ICL and the ML estimate of the scale parameter.
#'
logL.sparseBC <- function(x, results = NULL, Cs = NULL, Ds = NULL, mu = NULL, lambda = 0){
    if(attr(results, "class") == "sparseBC"){
        Cs <- results$Cs
        Ds <- results$Ds
        mu <- results$Mu
        lambda <- as.numeric(strsplit(as.character(results$cl), " ")[[5]])
        #cat(paste("taken lambda =",lambda,"\n"))
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
