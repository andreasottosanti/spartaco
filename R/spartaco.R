#' SpaRTaCo co-clustering
#'
#' This function returns the estimated model parameters and the co-clustering labels.
#'
#' @import SpatialExperiment
#' @export
#'
#' @param data either a `SpatialExperiment` object or a matrix containing the experiment;
#' @param coordinates if `is.matrix(data)`, it takes the matrix of spatial coordinates of dimension `ncol(data)` x 2.
#' @param assay if `class(data) == "SpatialExperiment"`, it takes either the name or the index of the assay to be used;
#' @param K the number of row clusters (only when `input.values == NULL`).
#' @param R the number of column clusters (only when `input.values == NULL`).
#' @param Delta.constr the constraint on the Delta matrix (default is 10; see **Details**).
#' @param max.iter the maximum number of iterations the estimation algorithm is run.
#' @param metropolis.iterations the number of iterations within each SE Step.
#' @param estimate.iterations the maximum number of iterations within each M Step.
#' @param prob.m the vector of probabilities assigned to the vector `1:length(prob.m)`, determining the number of columns the algorithm attempts to reallocate in a single SE Step.
#' @param input.values the starting points of the estimation process (see **Details**). If passed, the values of `K` and `R` are taken from it. The output of a previous model estimation can be passed here.
#' @param conv.criterion a list containing the parameters that define a converge criterion (see **Details**).
#' @param verbose if `TRUE`, it displays the on-going estimation process.
#' @param save.options a list for specifying the saving parameters (see **Details**).
#' @param seed set the interval seed of the function.
#'
#'
#' @return An object of class `spartaco` with the parameter estimates, the clustering labels, the log-likelihood value at each iteration and the data, the ICL, the data matrix and the coordinates matrix.

#'
#' @details `Delta.constr` gives the quantity \deqn{c = \tau_{kr} + \xi_{kr},} where \eqn{\tau_{kr}} and \eqn{\xi_{kr}} are the spatial variance and the nugget effect of block \eqn{(k,r)}.
#'
#' The algorithm can be initiated from a given set of starting values. To do so, `input.values` receives a list of the form
#' `list(mu, tau, xi, alpha, beta, phi, Cs, Ds)`, where:
#' - `mu`, `tau`, `xi`, `alpha` and `beta` are `K` x `R` matrices;
#' - `phi` is a vector of length `R`;
#' - `Cs` is a vector of length `nrow(data)` containing the row clustering labels;
#' - `Ds` is a vector of length `ncol(data)` containing the column clustering labels.
#'
#' If the algorithm is initiated from some starting values,  `K` and `R` are set automatically according to the input values.
#' If an object of class `spartaco` is passed to `input.values`, the estimation starts from the final estimate of the previous run (see **Examples**).
#'
#' If `conv.criterion == NULL`, the algorithm is stopped after `max.iter` itereations, otherwise it is stopped when the increment of the log-likelihood is smaller than a certain threshold `conv.criterion$epsilon` for `conv.criterion$iterations` times in a row.
#'
#' The function allows also to save the results even if the estimation is not completed. To do so, `save.options` receives a list of two parameters:
#' `after` gives the number of iterations after which the results are saved, `file.name` contains the path where the results are saved.
#'
# If `verbose == "full"`, the on-going estimation procedure is displayed. If `verbose == "progress"`, a dynamic progress bar will display the percentage of iterations completed.
# If `verbose == F`, then nothing is displayed in console.

# \deqn{p(x) = \frac{\lambda^x e^{-\lambda}}{x!}}{%p(x) = \lambda^x exp(-\lambda)/x!} for \eqn{x = 0, 1, 2, \ldots}

#' @references
#' Sottosanti, A. and Risso, D. (2021+) Co-clustering of Spatially Resolved Transcriptomic Data [(preprint)](https://arxiv.org/abs/2110.04872)

#' @examples
#' library(spartaco)
#'
#' # When x is a matrix:
#' x <- matrix(runif(n*p), n, p)
#' coordinates <- matrix(runif(2*p), p, 2)
#' output <- spartaco(data = x, coordinates = coordinates, K = K, R = R)
#'
#' # When x is an object of class "SpatialExperiment":
#' library(spatialLIBD)
#' ehub <- ExperimentHub::ExperimentHub()
#' if (!exists("spe")) spe <- fetch_data(type = "spe", eh = ehub)
#' output <- spartaco(data = spe, assay = 2, K = K, R = R)
#'
#' # To start the algorithm from the output of a previous run
#' output2 <- spartaco(data = x, coordinates = coordinates, input.val = output)

spartaco <- function(data,
                     coordinates = NULL,
                     assay = NULL,
                     K = NULL,
                     R = NULL,
                     Delta.constr = 10,
                     max.iter = 1000,
                     metropolis.iterations = 150,
                     estimate.iterations = 100,
                     prob.m = c(.7, .2, .1),
                     input.values = NULL,
                     conv.criterion = list(iterations = 10, epsilon = 1e-4),
                     verbose = T,
                     save.options = NULL,
                     seed = NULL
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
         prob.m = prob.m,
         conv.criterion = conv.criterion,
         input.values = input.values,
         save.options = save.options,
         verbose = verbose)
}
