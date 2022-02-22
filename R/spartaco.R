#' Spartaco co-clustering
#'
#' This function will return the estimated parameters and the co-cluster labels.
#'
#' @import SpatialExperiment
#' @export
#'
#' @param x the input data matrix.
#' @param coordinates the matrix of spatial coordinates of dimension ncol(x) X 2.
#' @param K the number of row clusters (when `input.values = NULL`).
#' @param R the number of column clusters (when `input.values = NULL`).
#' @param Delta.constr the constraint on the Delta matrix (default is 10; see Details).
#' @param max.iter the maximum number of iterations the estimation algorithm is run.
#' @param metropolis.iterations the number of iterations within each SE Step.
#' @param estimate.iterations the maximum number of iterations within each M Step.
#' @param sampling.m the method used for selecting *m*, the number of elements that the SE Step attempts to change in a single iteration (see Details).
#' @param prob.m the vector of probabilities used to sample *m* when `sampling.m = "standard"`.
#' @param input.values the starting points of the estimation process (see Details). If passed, the values of `K` and `R` are taken from it. The output of a previous model estimation can be passed here.
#' @param conv.criterion a list containing the parameters that define a converge criterion (see Details).
#' @param verbose if `TRUE`, it displays the on-going estimation process.
#' @param seed set the interval seed of the function.
#'
#' @details `Delta.constr` gives the quantity $\tau_{kr}$

#' @examples
#' library(SpatialExperiment)
#' example(SpatialExperiment)
#' # spartaco(se)
spartaco <- function(x,
                     coordinates,
                     K = NULL,
                     R = NULL,
                     Delta.constr = 10,
                     max.iter = 1000,
                     metropolis.iterations = 150,
                     estimate.iterations = 10,
                     sampling.m = "standard",
                     prob.m = c(.7, .2, .1),
                     input.values = NULL,
                     conv.criterion = list(iterations = 10, epsilon = 1e-4),
                     verbose = FALSE,
                     seed = NULL
                     ) {

    set.seed(seed = seed)
    Dist <- as.matrix(stats::dist(coordinates))

    if(!is.null(input.values)){
        K <- nrow(input.values$mu)
        R <- ncol(input.values$mu)
    }

    main(x = x,
         Dist = Dist,
         K = K, R = K,
         Delta.constr = Delta.constr,
         max.iter = max.iter,
         metropolis.iterations = metropolis.iterations,
         estimate.iterations = estimate.iterations,
         sampling.m = sampling.m,
         prob.m = prob.m,
         conv.criterion = conv.criterion,
         input.values = input.values,
         verbose = verbose)
}
