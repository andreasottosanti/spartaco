#' Multiple runs of SpaRTaCo
#'
#' This function returns the estimated model parameters and the co-clustering labels obtained after running SpaRTaCo multiple times in parallel.
#'
#' @import SpatialExperiment
#' @import parallel
#' @export
#'
#' @param data either a `SpatialExperiment` object or a matrix containing the experiment;
#' @param coordinates if `is.matrix(data)`, it takes the matrix of spatial coordinates of dimension `ncol(data)` x 2.
#' @param assay if `class(data) == "SpatialExperiment"`, it takes either the name or the index of the assay to be used;
#' @param K the number of row clusters (only when `input.values == NULL`);
#' @param R the number of column clusters (only when `input.values == NULL`);
#' @param nstart the number of parallel runs of the estimation algorithm;
#' @param Delta.constr the constraint on the Delta matrix (default is 10; see **Details**).
#' @param max.iter the maximum number of iterations the estimation algorithm is run.
#' @param metropolis.iterations the number of iterations within each SE Step.
#' @param estimate.iterations the maximum number of iterations within each M Step.
#' @param prob.m the vector of probabilities assigned to the vector `1:length(prob.m)`, determining the number of columns the algorithm attempts to reallocate in a single SE Step.
#' @param conv.criterion a list containing the parameters that define a converge criterion (see **Details**).
#' @param mc.cores the number of cores to be used.
#'
#' @return An object of class `spartaco` with the parameter estimates, the clustering labels, the log-likelihood value at each iteration and the data, the ICL, the data matrix and the coordinates matrix, and the clustering uncertainty.
#'
#' @details For details about the parameters, check [spartaco::spartaco()]

#'
#' @references
#' Sottosanti, A. and Risso, D. (2021+) Co-clustering of Spatially Resolved Transcriptomic Data [(preprint)](https://arxiv.org/abs/2110.04872)


spartaco_multirun <- function(data,
                     coordinates = NULL,
                     assay = NULL,
                     K = NULL,
                     R = NULL,
                     nstart = 5,
                     Delta.constr = 10,
                     max.iter = 1000,
                     metropolis.iterations = 150,
                     estimate.iterations = 100,
                     prob.m = c(.7, .2, .1),
                     conv.criterion = list(iterations = 10, epsilon = 1e-4),
                     mc.cores = 1
) {

    if(class(data) == "SpatialExperiment"){
        if(is.numeric(assay)) which.assay <- assay
        else which.assay <- which(names(data@assays@data) == assay)
        x <- as.matrix(data@assays@data[[which.assay]])
        row.names(x) <- rowData(data)$gene_name
        coordinates <- as.matrix(spatialCoords(data))
    } else {
        x <- data
    }

    results <- mclapply(1:nstart, function(l)
        spartaco(data = x, coordinates = coordinates, K = K, R = R, assay = NULL,
                 Delta.constr = Delta.constr, max.iter = max.iter,
                 metropolis.iterations = metropolis.iterations,
                 estimate.iterations = estimate.iterations,
                 prob.m = prob.m, conv.criterion = conv.criterion,
                 verbose = F, save.options = NULL), mc.cores = mc.cores)
    output <- CombineSpartaco(results)
    return(output)
}
