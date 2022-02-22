RowClustering <- function(x, Ds, Mu, Tau, Xi, Alpha, Beta, Phi, Uglob, Dglob){
    K <- ifelse(is.vector(Mu), 1, nrow(Mu))
    R <- length(Phi)
    ll <- matrix(0, nrow(x), K)
    goodR <- sort(unique(Ds))
    sapply(goodR, function(r){
        bl1 <- x[,Ds == r] %*% Uglob[[r]]
        bl2 <- matrix(1, nrow(x), sum(Ds == r)) %*% Uglob[[r]]
        sapply(1:K, function(k){
            alpha.i <- sum(Ds == r)/2 + Alpha[k,r]
            log.determ <- -.5 * sum(log(Tau[k,r]*Dglob[Ds == r] + Xi[k,r]))
            Block1 <- bl1 - Mu[k,r] * bl2
            beta.i <- (Block1 * Block1) %*% (1/(Tau[k,r]*Dglob[Ds == r] + Xi[k,r]))/2 + Beta[k,r]
            ll[,k] <<- ll[,k]+(
                -sum(Ds == r)/2*log(2*pi)+log.determ+
                    (Alpha[k,r]*log(Beta[k,r])-
                         lgamma(Alpha[k,r]))+
                    (lgamma(alpha.i)-alpha.i*log(beta.i))
            )
        })
    })
    allocation <- apply(ll,1,which.max)
    #stoch.allocation <- apply(ll,1,function(p){
    #    pp <- p-mean(p)
    #    sample(1:K, size = 1, prob = exp(pp))})
    return(list(allocation = allocation,
                #stoch.allocation = stoch.allocation,
                prob = ll))
}
