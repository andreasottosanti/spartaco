
GibbsAllocation <- function(x,
                                 Cs, Ds,
                                 Uglob, # list of matrices of length R, contains the matrices of eigenvec
                                 Dglob, # vector of length p, contains the eigenval
                                 Dist,
                                 Mu, Tau, Xi, Alpha, Beta, Phi,
                                 maxit = 10,
                                 min.obs = 3,
                                 prob.m = c(.7, .2, .1)){
    K <- ifelse(is.vector(Mu), 1, nrow(Mu))
    R <- length(Phi)
    if(is.vector(Mu)) Mu <- matrix(Mu, K, R)
    if(is.vector(Tau)) Tau <- matrix(Tau, K, R)
    if(is.vector(Xi)) Xi <- matrix(Xi, K, R)
    if(is.vector(Alpha)) Alpha <- matrix(Alpha, K, R)
    if(is.vector(Beta)) Beta <- matrix(Beta, K, R)
    goodK <- sort(unique(Cs))
    goodR <- sort(unique(Ds))
    logL.values <- matrix(0,K,R)
    sapply(goodR, function(r){
        sapply(goodK, function(k){
            logL.values[k,r] <<- logL.Cocluster(x = x[Cs == k, Ds == r], Mu = Mu[k,r], Tau = Tau[k,r], Xi = Xi[k,r], Alpha = Alpha[k,r],
                                                Beta = Beta[k,r], U = Uglob[[r]], d = Dglob[Ds == r])
        })
    })
    Ds.star <- Ds
    Uglob.star <- Uglob
    Dglob.star <- Dglob
    for(j in 1:ncol(x)){
        probabilities <- numeric(R)
        # move a single observation
        probabilities[Ds[j]] <- sum(logL.values)
        for(r in setdiff(1:R, Ds[j])){
            # initialization
            logL.values.star <- logL.values
            Ds.star <- Ds
            Uglob.star <- Uglob
            Dglob.star <- Dglob
            # computation
            Ds.star[j] <- r
            to.be.changed <- c(Ds[j], r)
            sapply(to.be.changed, function(r){
                eigK <- eigen(exp(-Dist[Ds.star == r, Ds.star == r]/Phi[r]))
                Uglob.star[[r]] <<- eigK$vec
                Dglob.star[Ds.star == r] <<- eigK$val
                sapply(goodK, function(k){
                    logL.values.star[k,r] <<- logL.Cocluster(x = x[Cs == k, Ds.star == r], Mu = Mu[k,r], Tau = Tau[k,r], Xi = Xi[k,r], Alpha = Alpha[k,r],
                                                             Beta = Beta[k,r], U = Uglob.star[[r]], d = Dglob.star[Ds.star == r])
                })
            })
            probabilities[r] <- sum(logL.values.star)
        }
        to.be.changed <- c(Ds[j], r)
        Ds[j] <- which.max(probabilities)
        sapply(to.be.changed, function(r){
            eigK <- eigen(exp(-Dist[Ds == r, Ds == r]/Phi[r]))
            Uglob[[r]] <<- eigK$vec
            Dglob[Ds == r] <<- eigK$val
            sapply(goodK, function(k){
                logL.values[k,r] <<- logL.Cocluster(x = x[Cs == k, Ds == r], Mu = Mu[k,r], Tau = Tau[k,r], Xi = Xi[k,r], Alpha = Alpha[k,r],
                                                         Beta = Beta[k,r], U = Uglob[[r]], d = Dglob[Ds == r])
            })
        })
    }

    results <- list(Ds = Ds, Uglob = Uglob, Dglob = Dglob, logL = sum(logL.values), logL.values = logL.values)
    return(results)
}



