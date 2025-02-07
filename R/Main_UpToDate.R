main <- function(x, coordinates,
                 K,
                 R,
                 Delta.constr = 10,
                 max.iter = 10^3,
                 constant.alpha.beta = F,
                 column.allocation.algorithm = "CEM",
                 metropolis.iterations = 150,
                 estimate.iterations = 10,
                 prob.m = c(.7, .2, .1),
                 conv.criterion = NULL,
                 input.values = NULL,
                 save.options = NULL,
                 verbose = F){
    Dist <- as.matrix(stats::dist(coordinates))
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
    #if(verbose == "progress"){
    #    progressr::handlers(global = T)
    #    P <- progressr::progressor(along = 1:max.iter)}
    while(T){
        if(i == max.iter) break
        i <- i + 1
        #if(verbose == "progress") P()
        if(verbose == T) cat(paste("---Iteration",i,"\n"))

        # ---M Step
        if(verbose == T) cat("M Step/")
        goodK <- sort(unique(cur.Cs))
        goodR <- sort(unique(cur.Ds))
        if(!constant.alpha.beta){
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
                                                                                                     maxit = estimate.iterations
                                                                                                     )
                    cur.mu[k,r] <<- estimation.parameters$mu
                    cur.tau[k,r] <<- estimation.parameters$tau
                    cur.xi[k,r] <<- estimation.parameters$xi
                    cur.alpha[k,r] <<- estimation.parameters$alpha
                    cur.beta[k,r] <<- estimation.parameters$beta
                })
            })
        } else {
            sapply(goodR, function(r){
                traceDelta_r <- Delta.constr * sum(cur.Ds == r)
                sapply(goodK, function(k){
                    estimation.parameters <- Estimate.Cocluster.Parameters.marginal.fixed.alpha.beta(x = x[cur.Cs == k, cur.Ds == r],
                                                                                                     traceDelta = traceDelta_r,
                                                                                                     U = Uglob[[r]],
                                                                                                     d = Dglob[cur.Ds == r],
                                                                                                     mu0 = cur.mu[k,r],
                                                                                                     alpha0 = cur.alpha[k,r],
                                                                                                     beta0 = cur.beta[k,r],
                                                                                                     tau0 = cur.tau[k,r],
                                                                                                     maxit = estimate.iterations
                    )
                    cur.mu[k,r] <<- estimation.parameters$mu
                    cur.tau[k,r] <<- estimation.parameters$tau
                    cur.xi[k,r] <<- estimation.parameters$xi
                })
            })
            alpha.beta <- Estimation.alpha.beta.fixed(alpha0 = unique(cur.alpha)[1,1],
                                                      beta0 = unique(cur.beta)[1,1],
                                                      x = x,
                                                      U = Uglob,
                                                      d = Dglob,
                                                      mu = cur.mu,
                                                      tau= cur.tau,
                                                      xi = cur.xi,
                                                      Cs = cur.Cs,
                                                      Ds = cur.Ds)
            cur.alpha <- matrix(alpha.beta[1], K, R)
            cur.beta <- matrix(alpha.beta[2], K, R)
        }

        # ---SE Step
        if(verbose == T) cat("SE Step/")
        if(column.allocation.algorithm == "SE")
            {cur.ds <- tryCatch({
                MetropolisAllocation(x = x, Uglob = Uglob, Dglob = Dglob,
                                               Cs = cur.Cs, Ds = cur.Ds, Dist = Dist, Mu = cur.mu, Tau = cur.tau, Xi = cur.xi, Alpha = cur.alpha, Beta = cur.beta, Phi = cur.phi,
                                               maxit = metropolis.iterations,
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
                })} else
                {
                    cur.ds <- GibbsAllocation(x = x, Cs = cur.Cs, Ds = cur.Ds, Uglob = Uglob, Dglob = Dglob, Dist = Dist, Mu = cur.mu, Tau = cur.tau, Xi = cur.xi, Alpha = cur.alpha, Beta = cur.beta, Phi = cur.phi)
                }
        cur.Ds <- cur.ds$Ds
        Uglob <- cur.ds$Uglob
        Dglob <- cur.ds$Dglob
        logL.values <- cur.ds$logL.values
        ll[i] <- sum(logL.values)

        # ---M Step
        if(verbose == T) cat("M Step/")
        goodK <- sort(unique(cur.Cs))
        goodR <- sort(unique(cur.Ds))
        sapply(goodR, function(r){
            traceDelta_r <- Delta.constr * sum(cur.Ds == r)
            sapply(goodK, function(k){
                if(!constant.alpha.beta){
                    estimation.parameters <- Estimate.Cocluster.Parameters.marginal.constraint.trace(x = x[cur.Cs == k, cur.Ds == r],
                                                                                                     traceDelta = traceDelta_r,
                                                                                                     U = Uglob[[r]],
                                                                                                     d = Dglob[cur.Ds == r],
                                                                                                     mu0 = cur.mu[k,r],
                                                                                                     alpha0 = cur.alpha[k,r],
                                                                                                     beta0 = cur.beta[k,r],
                                                                                                     tau0 = cur.tau[k,r],
                                                                                                     maxit = estimate.iterations
                    )
                    cur.alpha[k,r] <<- estimation.parameters$alpha
                    cur.beta[k,r] <<- estimation.parameters$beta
                } else {
                    estimation.parameters <- Estimate.Cocluster.Parameters.marginal.constraint.trace(x = x[cur.Cs == k, cur.Ds == r],
                                                                                                     traceDelta = traceDelta_r,
                                                                                                     U = Uglob[[r]],
                                                                                                     d = Dglob[cur.Ds == r],
                                                                                                     mu0 = cur.mu[k,r],
                                                                                                     alpha0 = cur.alpha[k,r],
                                                                                                     beta0 = cur.beta[k,r],
                                                                                                     tau0 = cur.tau[k,r],
                                                                                                     maxit = estimate.iterations
                    )
                }
                cur.mu[k,r] <<- estimation.parameters$mu
                cur.tau[k,r] <<- estimation.parameters$tau
                cur.xi[k,r] <<- estimation.parameters$xi
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
        if(constant.alpha.beta){
            alpha.beta <- Estimation.alpha.beta.fixed(alpha0 = unique(cur.alpha),
                                                      beta0 = unique(cur.beta),
                                                      x = x,
                                                      U = Uglob,
                                                      d = Dglob,
                                                      mu = cur.mu,
                                                      tau= cur.tau,
                                                      xi = cur.xi,
                                                      Cs = cur.Cs,
                                                      Ds = cur.Ds)
            cur.alpha <- matrix(alpha.beta[1], K, R)
            cur.beta <- matrix(alpha.beta[2], K, R)
        }

        # ---CE Step
        if(verbose == T) cat("CE Step/")
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

        if(verbose == T){
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

        # save the result in the given location
        if(!is.null(save.options)){
            if(i %% save.options$after == 0){
                ICL <- max(ll) - nrow(x)*K - ncol(x)*R - .5*(4*K*R+R)*log(nrow(x) * ncol(x))
                results <- list(
                    mu = best.mu,
                    tau = best.tau,
                    xi = best.xi,
                    alpha = best.alpha,
                    beta = best.beta,
                    phi = best.phi,
                    Cs = best.Cs,
                    Ds = best.Ds,
                    logL = ll[c(2:i)],
                    ICL = ICL,
                    x = x,
                    coordinates = coordinates
                )
                class(results) <- "spartaco"
                save(results, file = save.options$file.name)
            }
        }
    }

    ICL <- max(ll) - nrow(x)*log(K) - ncol(x)*log(R) - .5*(4*K*R+R)*log(nrow(x) * ncol(x))
    # save the result in the given location
    if(!is.null(save.options)){
        results <- list(
            mu = best.mu,
            tau = best.tau,
            xi = best.xi,
            alpha = best.alpha,
            beta = best.beta,
            phi = best.phi,
            Cs = best.Cs,
            Ds = best.Ds,
            logL = ll[c(2:i)],
            ICL = ICL,
            x = x,
            coordinates = coordinates)
        class(results) <- "spartaco"
        save(results, file = save.options$file.name)
    }
    results <- list(
        mu = best.mu,
        tau = best.tau,
        xi = best.xi,
        stn.ratio = best.tau/best.xi,
        alpha = best.alpha,
        beta = best.beta,
        phi = best.phi,
        Cs = best.Cs,
        Ds = best.Ds,
        logL = ll[c(2:i)],
        ICL = ICL,
        x = x,
        coordinates = coordinates
    )
    class(results) <- "spartaco"
    return(results)
}
