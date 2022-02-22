main <- function(x, Dist,
                 K,
                 R,
                 Delta.constr = 10,
                 max.iter = 10^3,
                 metropolis.iterations = 150,
                 estimate.iterations = 10,
                 sampling.m = "standard",
                 prob.m = c(.7, .2, .1),
                 conv.criterion = conv.criterion,
                 input.values = NULL,
                 verbose = F){
    if(is.null(input.values)){
        cur.Cs <- best.Ds <- sample(1:K, size = nrow(x), replace = T)
        cur.Ds <- best.Cs <- sample(1:R, size = ncol(x), replace = T)
        cur.phi <- best.phi <- runif(R, 1, 5)
        cur.mu <- best.mu <- matrix(runif(K*R, 1, 10), K, R)
        cur.tau <- best.tau <- matrix(runif(K * R, 1e-7, Delta.constr), K, R)
        cur.xi <- best.xi <- Delta.constr - best.tau
        cur.alpha <- best.alpha <- matrix(runif(K*R, 1, 3), K, R)
        cur.beta <- best.beta <- matrix(runif(K*R, 1, 3), K, R)} else {
            cur.Cs <- best.Cs <- input.values$Cs
            cur.Ds <- best.Ds <- input.values$Ds
            cur.phi <- best.phi <- input.values$phi
            cur.mu <- best.mu <- input.values$mu
            cur.tau <- best.tau <- input.values$tau
            cur.xi <- best.xi <- Delta.constr - best.tau
            cur.alpha <- best.alpha <- input.values$alpha
            cur.beta <- best.beta <- input.values$beta}
    Cs <- matrix(0, nrow(x), max.iter)
    Cs[,1] <- cur.Cs
    Ds <- matrix(0, ncol(x), max.iter)
    Ds[,1] <- cur.Ds

    Uglob <- list()
    Dglob <- numeric(ncol(x))
    sapply(1:R, function(r){
        eigK <- eigen(exp(-Dist[cur.Ds == r, cur.Ds == r]/cur.phi[r]))
        Uglob[[r]] <<- eigK$vec
        Dglob[cur.Ds == r] <<- eigK$val
    })
    ll <- rep(-1e+40, max.iter)
    logL.values <- matrix(0, K, R)
    i <- 1
    if(!is.null(conv.criterion)) counter.conv <- 0
    while(T){
        if(i == max.iter) break
        i <- i+1
        if(verbose) cat(paste("---Iteration",i,"\n"))

        # ---M Step
        if(verbose) cat("M Step/")
        goodK <- sort(unique(cur.Cs))
        goodR <- sort(unique(cur.Ds))
        sapply(goodR, function(r){
            traceDelta_r <- Delta.constr * sum(cur.Ds == r)
            sapply(goodK, function(k){
                estimation.parameters <- Estimate.Cocluster.Parameters.marginal.constraint.trace(x = x[cur.Cs == k, cur.Ds == r],
                                                                                                 traceDelta = traceDelta_r,
                                                                                                 U = Uglob[[r]],
                                                                                                 d = Dglob[cur.Ds == r],
                                                                                                 mu0 = cur.mu[k,r],
                                                                                                 alpha0 = cur.alpha[k,r],
                                                                                                 beta0 = cur.beta[k,r],
                                                                                                 tau0 = cur.tau[k,r],
                                                                                                 maxit = estimate.iterations)
                cur.mu[k,r] <<- estimation.parameters$mu
                cur.tau[k,r] <<- estimation.parameters$tau
                cur.xi[k,r] <<- estimation.parameters$xi
                cur.alpha[k,r] <<- estimation.parameters$alpha
                cur.beta[k,r] <<- estimation.parameters$beta
            })
            cur.phi[r] <<- updatePhi_r_marginal(x = x[,cur.Ds == r],
                                               Cs = cur.Cs,
                                               Dist = Dist[cur.Ds == r, cur.Ds == r],
                                               Mu = cur.mu[,r],
                                               Tau = cur.tau[,r],
                                               Xi = cur.xi[,r],
                                               Alpha = cur.alpha[,r],
                                               Beta = cur.beta[,r],
                                               phi.old = cur.phi[r])
            EigenK <- eigen(exp(-Dist[cur.Ds == r, cur.Ds == r]/cur.phi[r]))
            Uglob[[r]] <<- EigenK$vec
            Dglob[cur.Ds == r] <<- EigenK$val
        })

        # ---SE Step
        if(verbose) cat("SE Step/")
        cur.ds <- tryCatch({
            MetropolisAllocation(x = x, Uglob = Uglob, Dglob = Dglob,
                                           Cs = cur.Cs, Ds = cur.Ds, Dist = Dist, Mu = cur.mu, Tau = cur.tau, Xi = cur.xi, Alpha = cur.alpha, Beta = cur.beta, Phi = cur.phi,
                                           maxit = metropolis.iterations,
                                 sampling.m = sampling.m,
                                 rate.m = 1/(i-1)+.5,
                                 prob.m = prob.m,
                                 min.obs = 10)
            },
            error = function(cond){
                message(cond)
                llv <- matrix(0, K, R)
                sapply(goodR, function(r){
                    sapply(goodK, function(k){
                        llv[k,r] <- logL.Cocluster(x = x[cur.Cs == k, cur.Ds == r],
                                                            Mu = cur.mu[k,r],
                                                            Tau = cur.tau[k,r],
                                                            Xi = cur.xi[k,r],
                                                            Alpha = cur.alpha[k,r],
                                                            Beta = cur.beta[k,r],
                                                            U = Uglob[[r]],
                                                            d = Dglob[cur.Ds == r])
                    })
                })
                return(list(Ds = cur.Ds, Uglob = Uglob, Dglob = Dglob, logL.values = llv, accepted = 0))
            })
        cur.Ds <- cur.ds$Ds
        Uglob <- cur.ds$Uglob
        Dglob <- cur.ds$Dglob
        logL.values <- cur.ds$logL.values
        ll[i] <- sum(logL.values)

        # ---M Step
        if(verbose) cat("M Step/")
        goodK <- sort(unique(cur.Cs))
        goodR <- sort(unique(cur.Ds))
        sapply(goodR, function(r){
            traceDelta_r <- Delta.constr * sum(cur.Ds == r)
            sapply(goodK, function(k){
                estimation.parameters <- Estimate.Cocluster.Parameters.marginal.constraint.trace(x = x[cur.Cs == k, cur.Ds == r],
                                                                                                 traceDelta = traceDelta_r,
                                                                                                 U = Uglob[[r]],
                                                                                                 d = Dglob[cur.Ds == r],
                                                                                                 mu0 = cur.mu[k,r],
                                                                                                 alpha0 = cur.alpha[k,r],
                                                                                                 beta0 = cur.beta[k,r],
                                                                                                 tau0 = cur.tau[k,r],
                                                                                                 maxit = estimate.iterations)
                cur.mu[k,r] <<- estimation.parameters$mu
                cur.tau[k,r] <<- estimation.parameters$tau
                cur.xi[k,r] <<- estimation.parameters$xi
                cur.alpha[k,r] <<- estimation.parameters$alpha
                cur.beta[k,r] <<- estimation.parameters$beta
            })
            cur.phi[r] <<- updatePhi_r_marginal(x = x[,cur.Ds == r],
                                                Cs = cur.Cs,
                                                Dist = Dist[cur.Ds == r, cur.Ds == r],
                                                Mu = cur.mu[,r],
                                                Tau = cur.tau[,r],
                                                Xi = cur.xi[,r],
                                                Alpha = cur.alpha[,r],
                                                Beta = cur.beta[,r],
                                                phi.old = cur.phi[r])
            EigenK <- eigen(exp(-Dist[cur.Ds == r, cur.Ds == r]/cur.phi[r]))
            Uglob[[r]] <<- EigenK$vec
            Dglob[cur.Ds == r] <<- EigenK$val
        })

        # ---CE Step
        if(verbose) cat("CE Step/")
        cur.cs <- RowClustering(x = x, Ds = cur.Ds, Mu = cur.mu, Tau = cur.tau, Xi = cur.xi, Alpha = cur.alpha, Beta = cur.beta, Phi = cur.phi, Uglob = Uglob, Dglob = Dglob)
        cur.Cs <- cur.cs$allocation
        goodK <- sort(unique(cur.Cs))
        goodR <- sort(unique(cur.Ds))
        sapply(goodR, function(r){
            sapply(goodK, function(k){
                logL.values[k,r] <<- logL.Cocluster(x = x[cur.Cs == k, cur.Ds == r],
                                                           Mu = cur.mu[k,r],
                                                           Tau = cur.tau[k,r],
                                                           Xi = cur.xi[k,r],
                                                           Alpha = cur.alpha[k,r],
                                                           Beta = cur.beta[k,r],
                                                           U = Uglob[[r]],
                                                           d = Dglob[cur.Ds == r])
            })
        })
        ll[i] <- sum(logL.values)
        Ds[,i] <- cur.Ds
        Cs[,i] <- cur.Cs

        if(ll[i] == max(ll)){
            best.phi <- cur.phi
            best.mu <- cur.mu
            best.tau <- cur.tau
            best.xi <- cur.xi
            best.alpha <- cur.alpha
            best.beta <- cur.beta
            best.Cs <- cur.Cs
            best.Ds <- cur.Ds
        }

        if(verbose){
            cat(paste("diff(logL) =",round(diff(ll)[i-1],5),"\n"))
            cat(paste("Row clusters size =", paste(table(cur.Cs), collapse = ", "),"\n"))
            cat(paste("Column clusters size =", paste(table(cur.Ds), collapse = ", "),"\n"))
        }

        if(!is.null(conv.criterion)){
            if(ll[i] >= ll[i-1] & (ll[i] - ll[i-1] < conv.criterion$epsilon)){
                counter.conv <- counter.conv + 1
                if(counter.conv == conv.criterion$iterations){
                    cat("Converged\n")
                    break}
            } else {
                counter.conv <- 0
            }
        }
    }
    return(list(
        phi = best.phi,
        mu = best.mu,
        tau = best.tau,
        xi = best.xi,
        alpha = best.alpha,
        beta = best.beta,
        Cs = best.Cs,
        Ds = best.Ds,
        CS = Cs,
        DS = Ds,
        logL = ll,
        x = x
    ))
}
