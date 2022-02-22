logL.Cocluster <- function(x, Mu, Tau, Xi, Alpha, Beta, U, d){
    if(is.vector(x)) x <- matrix(x, nrow = 1)
    Block1 <- (x - Mu) %*% U
    invD <- 1/(Tau*d + Xi)
    alpha.post.i <- ncol(x)/2 + Alpha
    beta.post.i <- (Block1 * Block1) %*% (invD)/2 + Beta
    -nrow(x)*ncol(x)/2*log(2*pi) +
        nrow(x)/2*sum(log(invD))+
        nrow(x)*(Alpha*log(Beta)-lgamma(Alpha))+
        sum(lgamma(alpha.post.i)-alpha.post.i*log(beta.post.i))
}

