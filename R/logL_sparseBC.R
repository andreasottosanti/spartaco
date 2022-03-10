#' Log-likelihood and ICL of the Sparse Biclustering model
#'
#' This function computes the log-likelihood and the ICL of the Sparse Biclustering model of Tan & Witten (2014).
#'
#' @export
#'
#' @param x the data matrix;
#' @param Cs the row clustering labels;
#' @param Ds the column clustering labels;
#' @param mu the `K`x`R` matrix containing the co-cluster centroids;
#' @param lambda the penalization parameter;
#' @param results it recieves an object of class `sparseBC` as input;
#' @param ICL logical; if true, it returns the integrated complete log-likelihood of the model.
#'
#' @return If `ICL = F`, it returns the log-likelihood of the model, otherwise it returns the ICL.
#'
logL.sparseBC <- function(x, Cs = NULL, Ds = NULL, mu = NULL, lambda = 0, results = NULL, ICL = F){
    if(attr(results, "class") == "sparseBC"){
        Cs <- results$Cs
        Ds <- results$Ds
        mu <- results$Mu
    }
    nk <- table(Cs)
    pr <- table(ds)
    K <- length(nk)
    R <- length(pr)
    sigma2 <- sum(sapply(1:K, function(k){
        sapply(1:R, function(r){
            sum((x[Cs == k, Ds == r]-mu[k,r])^2)
        })
    }))/sum(nk %*% t(pr))
    logL <- matrix(0, K, R)
    sapply(1:K, function(k){
        sapply(1:R, function(r){
            logL[k,r] <<- dmnv(x[Cs == k, Ds == r], mu = mu[k,r], sigma = diag(sigma2, nrow = nk[k]), delta = diag(1, nrow = pr[r]))
        })
    })
    if(ICL) return(2*sum(logL)-lambda*sum(abs(mu))-nrow(x)*log(K)-ncol(x)*log(R)-(K*R+1)/2*prod(dim(X))) else
        return(sum(logL) - lambda * sum(abs(mu)))
}
