#' PCA K-means model selection
#'
#' This function returns the estimated number of row and column clusters by first applying a PCA on the data,
#' and then determining the best K-means model based on a change point model fitted on the total within sum of squares.
#'
#' @import SpatialExperiment
#' @export
#'
#' @param x the data matrix;
#' @param K.range the range of row clusters to evaluate;
#' @param R.range the range of column clusters to evaluate.
#'
#' @return a vector of length 2 containing the values of `K` and `R` selected.
#'
PCA.Kmeans.KR <- function(x, K.range, R.range){
    change.point.model <- function(y, x, nstart = 10, plot.fit = F){
    values <- matrix(0,nstart,5)
    for(n in 1:nstart){
        f <- function(x, param){
            x0 <- param[1]
            beta0 <- param[2]
            beta1 <- param[3]
            beta2 <- param[4]
            beta0 + (beta1*(x-x0))*I(x-x0 < 0) + beta2*I(x-x0 >= 0)
        }

        least.square.f <- function(x, y, x0, param){
            eta <- f(x = x, param = c(x0, param))
            sum((y-eta)^2)
        }

        min.value <- matrix(0, length(x), 4)
        for(k in 1:length(x)){
            opt <- optim(runif(3,0,10), fn = least.square.f, x = 1:length(x), y = y, x0 = x[k])
            min.value[k,1] <- opt$val
            min.value[k,-1] <- opt$par
        }

        values[n,] <- c(x[which.min(min.value[,1])], min.value[which.min(min.value[,1]),])
    }
    if(plot.fit){
        plot(x, y, type = "b")
        curve(f(x, param = values[which.min(values[,2]),-2]), col = 2, add = T)}

    return(values[which.min(values[,2]),1])
    }

    pc.row <- prcomp(x, center = T, retx = T)
    pc.col <- prcomp(t(x), center = T, retx = T)

    wssq.rows <- numeric(length(K.range))
    for(k in 1:length(K.range)){
        print(k)
        km.rows.k <- kmeans(x = pc.row$x, centers = K.range[k], nstart = 20)
        wssq.rows[k] <- km.rows.k$tot.withinss
    }
    K.sel <- change.point.model(y = wssq.rows, x = K.range, nstart = 20, plot.fit = T)

    wssq.cols <- numeric(length(R.range))
    for(r in 1:length(R.range)){
        print(r)
        km.cols.r <- kmeans(x = pc.col$x[,1:2], centers = R.range[r], nstart = 20)
        wssq.cols[r] <- km.cols.r$tot.withinss
    }
    R.sel <- change.point.model(y = wssq.cols, x = R.range, nstart = 20, plot.fit = T)
    return(c(K.values, R.values))
}
