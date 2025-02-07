change.point.model <- function(y, x, nstart = 10, continuity = F, plot.fit = F, plot.title = NULL){
    # internal functions
    if(continuity){
        f <- function(x, param){
            x0 <- param[1]
            beta0 <- param[2]
            beta1 <- param[3]
            beta0 + (beta1*(x-x0))*I(x-x0 < 0)}
    } else {
        f <- function(x, param){
            x0 <- param[1]
            beta0 <- param[2]
            beta1 <- param[3]
            beta2 <- param[4]
            beta0 + (beta1*(x-x0))*I(x-x0 < 0) + beta2*I(x-x0 >= 0)}
    }
    least.square.f <- function(x, y, x0, param){
        eta <- f(x = x, param = c(x0, param))
        sum((y-eta)^2)
    }

    # main of the function
    values <- matrix(0, nstart, 5)
    for(n in 1:nstart){
        min.value <- matrix(0, length(x), 4)
        for(k in 1:length(x)){
            opt <- optim(c(mean(y), runif(1,0,10), mean(y)), fn = least.square.f, x = 1:length(x), y = y, x0 = x[k])
            min.value[k,1] <- opt$val
            min.value[k,-1] <- opt$par
        }

        values[n,] <- c(x[which.min(min.value[,1])], min.value[which.min(min.value[,1]),])
    }
    if(plot.fit){
        plot(x, y, type = "b", main = plot.title, xlab = "N. of clusters", ylab = "within sum of squares")
        curve(f(x, param = values[which.min(values[,2]),-2]), col = 2, add = T)
        abline(v = values[which.min(values[,2]),1], lty = 2, col = 4)
    }

    return(values[which.min(values[,2]),1])
}
