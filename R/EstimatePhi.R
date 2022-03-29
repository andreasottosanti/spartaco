updatePhi_r_marginal <- function(x, Cs, Dist, Mu, Tau, Xi, Alpha, Beta, phi.old = 1, lambda.phi = NULL){
    if(is.null(lambda.phi)) lambda.phi <- 0
    goodK <- sort(unique(Cs))
    routine.phi <- optim(par = phi.old, fn = function(phi){
        if(phi <= 0) return(-Inf)
        eig <- eigen(exp(-Dist/phi))
        val <- sapply(goodK, function(k){
            block1 <- (x[Cs == k,] - Mu[k]) %*% eig$vec
            G.mat <- block1 * block1
            alpha.i <- ncol(x)/2 + Alpha[k]
            beta.i <- G.mat %*% (1/(Tau[k]*eig$val + Xi[k]))/2 + Beta[k]
            return(-sum(Cs == k)/2*sum(log(Tau[k]*eig$val + Xi[k])) - sum(alpha.i * log(beta.i)))
        })
        -sum(val - lambda.phi * phi^2)
    })
    if(routine.phi$conv != 0) stop("Converge error in Phi!")
    return(routine.phi$par)
}
