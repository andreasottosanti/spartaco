Estimate.Cocluster.Parameters.marginal.constraint.trace <- function(x,
                                                                    U,
                                                                    d,
                                                                    mu0,
                                                                    alpha0,
                                                                    beta0,
                                                                    tau0,
                                                                    traceDelta,
                                                                    lambda.tau = NULL,
                                                                    maxit = 200,
                                                                    threshold = 1e-4){
    if(is.null(lambda.tau)) lambda.tau <- 0
    n <- nrow(x)
    p <- ncol(x)
    Mu <- cur.mu <- mu0
    Tau <- cur.tau <- tau0
    cur.xi <- traceDelta/p - tau0
    Alpha <- cur.alpha <- alpha0
    Beta <- cur.beta <- beta0
    logL <- tryCatch(
        {
            logL.Cocluster(x, Mu, Tau, traceDelta/p-Tau, Alpha, Beta, U, d)
        },
        error = function(cond) {
            return(-1e-40)
        })
    converged <- F
    bl1 <- x %*% U
    bl2 <- matrix(1, n, p) %*% U
    cur.xi <- traceDelta/p - cur.tau
    for(i in 2:maxit){
        # --update mu
        A.mat <- bl1 %*% diag(1/(cur.tau * d + cur.xi)) %*% t(bl1)
        B.mat <- bl1 %*% diag(1/(cur.tau * d + cur.xi)) %*% t(bl2)
        C.mat <- bl2 %*% diag(1/(cur.tau * d + cur.xi)) %*% t(bl2)
        routine.mu <- optim(par = cur.mu, fn = function(mu){
            sum(
                log(
                    (diag(A.mat) - 2*mu*diag(B.mat) + mu^2*diag(C.mat))/2 + cur.beta
                )
            )*(p/2 + cur.alpha)
        })
        if(routine.mu$convergence != 0){
            stop("Convergence error in mu!")
        }
        cur.mu <- routine.mu$par

        # --update alpha
        quadratic <- (diag(A.mat) - 2*cur.mu*diag(B.mat) + cur.mu^2*diag(C.mat))/2
        routine.alpha <- optim(cur.alpha, function(a){
            if(a <= 0) return(-Inf)
            -(
                a*(n*log(cur.beta) - sum(log(quadratic + cur.beta)))-n*(lgamma(a)-lgamma(p/2+a))
            )})
        if(routine.alpha$convergence != 0){
            stop("Convergence error in alpha!")
        }
        cur.alpha <- routine.alpha$par

        # --update beta
        routine.beta <- optim(cur.beta, function(b){
            if(b <= 0) return(-Inf)
            -(n*cur.alpha*log(b)-(p/2+cur.alpha)*sum(log(quadratic+b)))})
        if(routine.beta$convergence != 0){
            stop("Convergence error in beta!")
        }
        cur.beta <- routine.beta$par

        # --update tau and xi
        Block1 <- bl1 - cur.mu * bl2
        starting.tau <- ifelse(cur.tau < traceDelta/p, cur.tau, runif(1, 1e-7, traceDelta/p))
        G.mat <- Block1 * Block1
        routine.tau <- optim(starting.tau,
                             fn = function(taup){
                                 if(taup <= 0) return(-Inf)
                                 xip <- traceDelta/p - taup
                                 if(xip <= 0) return(-Inf)
                                 -(
                                     -n/2*sum(log(taup * d + xip)) -
                                         (p/2+cur.alpha) * sum(log(G.mat %*% (1/(taup * d + xip))/2 + cur.beta))
                                     - lambda.tau * abs(taup)
                                    )
                             }, control = list(maxit = 1000))
        if(routine.tau$convergence != 0){
            stop("Convergence error in tau!")
        }
        cur.tau <- routine.tau$par
        cur.xi <- traceDelta/p - cur.tau

        Mu[i] <- cur.mu
        Alpha[i] <- cur.alpha
        Beta[i] <- cur.beta
        Tau[i] <- cur.tau

        logL[i] <- logL.Cocluster(x, Mu[i], Tau[i], traceDelta/p-Tau[i], Alpha[i], Beta[i], U, d)
        #if(round(logL[i] - logL[i-1], 2) < 0) stop(cat("Decreasing loglikelihood within the M Step:",logL[i-1],"and",logL[i]))
        if((logL[i] - logL[i-1]) < threshold){
            converged <- T
            break}
    }
    return(list(mu = Mu[i],
                alpha = Alpha[i],
                beta = Beta[i],
                tau = Tau[i],
                xi = cur.xi,
                logL = logL[i]))
}


