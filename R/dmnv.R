dmvn <- function(x, mu, sigma, delta, log.res = T){
    n <- nrow(x)
    p <- ncol(x)
    logL <- -n*p/2*log(2*pi)-
        p/2*determinant(sigma, logarithm = T)-
        n/2*determinant(delta, logarithm = T)-
        .5*sum(diag(solve(sigma) %*% (x-mu) %*% solve(delta) %*% t(x-mu)))
    if(log.res) return(logL) else return(exp(logL))
}
