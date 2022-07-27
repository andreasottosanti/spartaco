MetropolisAllocation2 <- function(x, Cs, Ds,
                                 Uglob, # list of matrices of length R, contains the matrices of eigenvec
                                 Dglob, # vector of length p, contains the eigenval
                                 Dist,
                                 Mu, Tau, Xi, Alpha, Beta, Phi, maxit = 10, min.obs = 3,
                                 rate.m = NULL,
                                 prob.m = c(.7, .2, .1)){
    K <- ifelse(is.vector(Mu), 1, nrow(Mu))
    R <- length(Phi)
    if(is.vector(Mu)) Mu <- matrix(Mu, K, R)
    if(is.vector(Tau)) Tau <- matrix(Tau, K, R)
    if(is.vector(Xi)) Xi <- matrix(Xi, K, R)
    if(is.vector(Alpha)) Alpha <- matrix(Alpha, K, R)
    if(is.vector(Beta)) Beta <- matrix(Beta, K, R)
    D <- as.vector(table(Ds))
    Ds.init <- Ds
    accepted <- F
    accepted <- m.list <- numeric(maxit)
    logL.values <- matrix(0,K,R)
    goodK <- sort(unique(Cs))
    goodR <- sort(unique(Ds))
    Delta.inv <-
    sapply(goodR, function(r){
        sapply(goodK, function(k){
            logL.values[k,r] <<- logL.Cocluster(x = x[Cs == k, Ds == r], Mu = Mu[k,r], Tau = Tau[k,r], Xi = Xi[k,r], Alpha = Alpha[k,r],
                                                Beta = Beta[k,r], U = Uglob[[r]], d = Dglob[Ds == r])
        })
    })
    logL.den <- sum(logL.values)
    for(iter in 1:maxit){
        m <- m.list[iter] <- sample(1:length(prob.m),1, prob = prob.m)
        move <- sample(c("M1","M2"), size = 1)
        if(move == "M1"){
            gr.start <- sample(1:length(D), 1)
            gr.end <- sample(setdiff(1:length(D), gr.start), 1)
            j <- sample(1:D[gr.start], m, replace = F)
            Ds.star <- Ds
            Ds.star[which(Ds == gr.start)][j] <- gr.end
            to.be.changed <- unique(c(gr.start, gr.end))
            logL.values.star <- logL.values
            Uglob.star <- Uglob
            Dglob.star <- Dglob
            sapply(to.be.changed, function(r){
                eigK <- eigen(exp(-Dist[Ds.star == r, Ds.star == r]/Phi[r]))
                Uglob.star[[r]] <<- eigK$vec
                Dglob.star[Ds.star == r] <<- eigK$val
                sapply(goodK, function(k){
                    logL.values.star[k,r] <<- logL.Cocluster(x = x[Cs == k, Ds.star == r], Mu = Mu[k,r], Tau = Tau[k,r], Xi = Xi[k,r], Alpha = Alpha[k,r],
                                                             Beta = Beta[k,r], U = Uglob.star[[r]], d = Dglob.star[Ds.star == r])
                })
            })
            log.proposal.num <- sum(log(D[gr.start]-0:(m-1)))
            log.proposal.den <- sum(log(D[gr.end]+1:m))
            logL.num <- sum(logL.values.star)
            A <- exp(logL.num + log.proposal.num - logL.den - log.proposal.den)*all(as.vector(table(Ds.star)) >= min.obs)
            if(runif(1) <= A){
                Ds <- Ds.star
                D <- as.vector(table(Ds))
                logL.values <- logL.values.star
                logL.den <- logL.num
                Uglob <- Uglob.star
                Dglob <- Dglob.star
                accepted[iter] <- 1}
        }
        if(move == "M2"){
            gr.start <- sample(1:R, m, replace = T)
            gr.end <- sapply(1:m, function(k) sample(setdiff(1:R, gr.start[k]), 1))
            q1r <- sapply(1:R, function(r) sum(gr.start == r))
            q2r <- sapply(1:R, function(r) sum(gr.end == r))
            j <- sapply(1:R, function(r) ifelse(q1r[r] != 0, return(sample(1:D[r], q1r[r], replace = F)), 0))
            Ds.star <- Ds
            for(r in which(q1r != 0)) Ds.star[which(Ds == r)][j[[r]]] <- gr.end[gr.start == r]
            to.be.changed <- unique(c(gr.start, gr.end))
            logL.values.star <- logL.values
            Uglob.star <- Uglob
            Dglob.star <- Dglob
            sapply(to.be.changed, function(r){
                eigK <- eigen(exp(-Dist[Ds.star == r, Ds.star == r]/Phi[r]))
                Uglob.star[[r]] <<- eigK$vec
                Dglob.star[Ds.star == r] <<- eigK$val
                sapply(goodK, function(k){
                    logL.values.star[k,r] <<- logL.Cocluster(x = x[Cs == k, Ds.star == r], Mu = Mu[k,r], Tau = Tau[k,r], Xi = Xi[k,r], Alpha = Alpha[k,r],
                                                             Beta = Beta[k,r], U = Uglob.star[[r]], d = Dglob.star[Ds.star == r])
                })
            })
            log.proposal.num <- sum(log(factorial(q2r)))+sum(sapply(which(q1r != 0), function(r) sum(log(D[r]+0:(q1r[r]-1)))))
            log.proposal.den <- sum(log(factorial(q1r)))+sum(sapply(which(q2r != 0), function(r) sum(log(D[r]-q1r[r]+1:q2r[r]))))
            logL.num <- sum(logL.values.star)
            A <- exp(logL.num + log.proposal.num - logL.den - log.proposal.den)*all(as.vector(table(Ds.star)) >= min.obs)
            if(runif(1) <= A){
                Ds <- Ds.star
                D <- as.vector(table(Ds))
                logL.values <- logL.values.star
                logL.den <- logL.num
                Uglob <- Uglob.star
                Dglob <- Dglob.star
                accepted[iter] <- 1}
        }
    }
    results <- list(Ds = Ds, Uglob = Uglob, Dglob = Dglob, logL = logL.den, logL.values = logL.values,
                    accepted = sum(Ds != Ds.init), m.list = m.list)
    return(results)
}
