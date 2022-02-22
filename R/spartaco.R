#' Spartaco co-clustering
#'
#' This function will return the estimated parameters and the co-cluster labels.
#'
#' @import SpatialExperiment
#' @export
#'
#' @examples
#' library(SpatialExperiment)
#' example(SpatialExperiment)
#' # spartaco(se)
spartaco <- function(x,
                     coordinates,
                     K = NULL,
                     R = NULL,
                     traceRatio = 10,
                     max.iter = 1000,
                     metropolis.iterations = 150,
                     estimate.iterations = 10,
                     sampling.m = "standard",
                     prob.m = c(.7, .2, .1),
                     input.values = NULL,
                     verbose = FALSE,
                     conv.criterion = list(iterations = 10, epsilon = 1e-4),
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
         traceRatio = traceRatio,
         max.iter = max.iter,
         metropolis.iterations = metropolis.iterations,
         estimate.iterations = estimate.iterations,
         sampling.m = sampling.m,
         prob.m = prob.m,
         conv.criterion = conv.criterion,
         input.values = input.values,
         verbose = verbose)
}
