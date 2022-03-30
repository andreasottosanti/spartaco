#' PCA K-means model selection
#'
#' This function returns the estimated number of row and column clusters by first applying a PCA on the data,
#' and then determining the best K-means model based on a change point model fitted on the total within sum of squares.
#'
# #' @import SpatialExperiment
#' @export
#'
#' @param x the data matrix;
#' @param K.max the maximum number of row clusters to evaluate;
#' @param R.max the maximum number of column clusters to evaluate.
#' @param continuity if `TRUE`, the change point linear model is fitted under the constraint that the function is continuous in the change point.
#'
#' @return A list containing the values of `K` and `R` selected and the clustering labels with the
#' selected dimensions.
#'
PCA.Kmeans.KR <- function(x, K.max, R.max, continuity = F){
    K.range <- 1:K.max
    R.range <- 1:R.max
    pc.row <- prcomp(x, center = T, retx = T)$x[,1:2]
    pc.col <- prcomp(t(x), center = T, retx = T)$x[,1:2]

    wssq.rows <- numeric(length(K.range))
    for(k in 1:length(K.range)){
        km.rows.k <- kmeans(x = pc.row, centers = K.range[k], nstart = 20, iter.max = 100)
        wssq.rows[k] <- km.rows.k$tot.withinss
    }
    K.sel <- change.point.model(y = wssq.rows, x = K.range, nstart = 20, continuity = continuity, plot.fit = T, plot.title = "N. of row clusters")
    row.cluster <- kmeans(x = pc.row, centers = K.sel, nstart = 20, iter.max = 100)$cluster

    wssq.cols <- numeric(length(R.range))
    for(r in 1:length(R.range)){
        km.cols.r <- kmeans(x = pc.col, centers = R.range[r], nstart = 20, iter.max = 100)
        wssq.cols[r] <- km.cols.r$tot.withinss
    }
    R.sel <- change.point.model(y = wssq.cols, x = R.range, nstart = 20, continuity = continuity, plot.fit = T, plot.title = "N. of column clusters")
    col.cluster <- kmeans(x = pc.col, centers = R.sel, nstart = 20, iter.max = 100)$cluster

    return(list(KR = c(K.sel, R.sel),
                Cs = row.cluster,
                Ds = col.cluster))
}
