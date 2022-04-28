#' SpaRTaCo residuals
#'
#' This function the residuals of a SpaRTaCo model.
#'
#' @export
#'
#' @param x a `spartaco` object.
#'
#' @return An object of class `spartaco` with the same values of the input object, expect `$x`, which now contains the matrix of the residuals.

residuals.spartaco <- function(x){
    if(class(x) != "spartaco") stop("x is not a spartaco object")
    K <- nrow(x$mu)
    R <- ncol(x$mu)
    res <- x
    sapply(1:R, function(r){
        sapply(1:K, function(k){
            Deltakr <- x$tau[k,r]*exp(-as.matrix(dist(x$coordinates[x$Ds == r,]))/x$phi[r])+diag(x$xi[k,r], nrow =  sum(x$Ds == r))
            Fkr <- t(chol(Deltakr))
            res$x[x$Cs == k, x$Ds == r] <<- (x$x[x$Cs == k, x$Ds == r]-x$mu[k,r]) %*% solve(Fkr)
        })
    })
    return(res)
}
