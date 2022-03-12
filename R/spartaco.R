#' SpaRTaCo co-clustering
#'
#' This function returns the estimated model parameters and the co-clustering labels.
#'
#' @import SpatialExperiment
#' @export
#'
#' @param x the input data matrix.
#' @param coordinates the matrix of spatial coordinates of dimension ncol(x) x 2.
#' @param K the number of row clusters (when `input.values = NULL`).
#' @param R the number of column clusters (when `input.values = NULL`).
#' @param Delta.constr the constraint on the Delta matrix (default is 10; see **Details**).
#' @param max.iter the maximum number of iterations the estimation algorithm is run.
#' @param metropolis.iterations the number of iterations within each SE Step.
#' @param estimate.iterations the maximum number of iterations within each M Step.
#' @param sampling.m the method used for selecting *m*, the number of elements that the SE Step attempts to change in a single iteration (see **Details**).
#' @param prob.m the vector of probabilities used to sample *m* when `sampling.m = "standard"`.
#' @param input.values the starting points of the estimation process (see **Details**). If passed, the values of `K` and `R` are taken from it. The output of a previous model estimation can be passed here.
#' @param conv.criterion a list containing the parameters that define a converge criterion (see **Details**).
#' @param verbose if `TRUE`, it displays the on-going estimation process.
#' @param save.options a list for specifying the saving parameters (see **Details**).
#' @param seed set the interval seed of the function.
#'
#'
#' @return The function returns a list with the parameter estimates, the clustering labels, the log-likelihood value at each iteration and the data, the ICL index, the data matrix and the coordinates matrix.

#'
#' @details `Delta.constr` gives the quantity \deqn{c = \tau_kr + \xi_kr,} where \eqn{\tau_kr} and \eqn{\xi_kr} are the spatial variance and the nugget effect of block \eqn{(k,r)}.
#'
#' The are two methods for selecting *m*, the number of items the SE Step attempts to change in a single iteration.
#' If `sampling.m = "standard"`, then `m = sample(1:length(prob.m), 1, F, prob.m)`.
#' If `sampling.m = "adaptive"`, then `m = 1 + rpois(1, 1/t+1/2)`, where `t` denotes the t-th iteration of the estimation algorithm.
#'
#' The algorithm can be initiated from a given set of starting values. To do so, `input.values` receives a list of the form
#' `list(mu, tau, xi, alpha, beta, phi, Cs, Ds)`, where:
#' - `mu`, `tau`, `xi`, `alpha` and `beta` are `K` x `R` matrices;
#' - `phi` is a vector of length `R`;
#' - `Cs` is a vector of length `nrow(x)` containing the row clustering labels;
#' - `Ds` is a vector of length `ncol(x)` containing the column clustering labels.
#'
#' If the algorithm is initiated from some starting values, it is no longer necessary to set `K` and `R`, as they are automatically taken from here.
#' By passing to `input.values` an output of the `spartaco` function, the algorithm starts the estimation from the arrival points of the previous run (see **Examples**).
#'
#' If `conv.criterion = NULL`, the algorithm is stopped after `max.iter` itereations, otherwise it is stopped when the increment of the log-likelihood is smaller than a certain threshold `conv.criterion$epsilon` for `conv.criterion$iterations` times in a row.
#'
#' The function allows also to save the results even if the estimation is not completed. To do so, `save.options` received a list of two parameters:
#' `after` gives the number of iterations after which the results are saved, `file.name` contains the path where the results are saved.

# \deqn{p(x) = \frac{\lambda^x e^{-\lambda}}{x!}}{%p(x) = \lambda^x exp(-\lambda)/x!} for \eqn{x = 0, 1, 2, \ldots}

#' @references
#' Sottosanti, A. and Risso, D. (2021) Co-clustering of Spatially Resolved Transcriptomic Data [(preprint)](https://arxiv.org/abs/2110.04872)

#' @examples
#' library(spartaco)
#' x <- matrix(runif(n*p), n, p)
#' coordinates <- matrix(runif(2*p), p, 2)
#' output <- spartaco(x = x, coordinates = coordinates, K = K, R = R)
#'
#' # To start the algorithm from the output of a previous run
#' output2 <- spartaco(x = x, coordinates = coordinates, input.val = output)
spartaco <- function(x,
                     coordinates,
                     K = NULL,
                     R = NULL,
                     Delta.constr = 10,
                     max.iter = 1000,
                     metropolis.iterations = 150,
                     estimate.iterations = 100,
                     sampling.m = "standard",
                     prob.m = c(.7, .2, .1),
                     input.values = NULL,
                     conv.criterion = list(iterations = 10, epsilon = 1e-4),
                     verbose = TRUE,
                     save.options = NULL,
                     seed = NULL
                     ) {

    set.seed(seed = seed)

    if(!is.null(input.values)){
        K <- nrow(input.values$mu)
        R <- ncol(input.values$mu)
    }

    main(x = x,
         coordinates = coordinates,
         K = K, R = R,
         Delta.constr = Delta.constr,
         max.iter = max.iter,
         metropolis.iterations = metropolis.iterations,
         estimate.iterations = estimate.iterations,
         sampling.m = sampling.m,
         prob.m = prob.m,
         conv.criterion = conv.criterion,
         input.values = input.values,
         save.options = save.options,
         verbose = verbose)
}
