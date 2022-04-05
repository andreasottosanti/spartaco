#' SpaRTaCo ICL
#'
#' This function returns the ICL of a SpaRTaCo model.
#'
# #' @import SpatialExperiment
#' @export
#'
#' @param x a `spartaco` object.
#'
#'@return The function returns the ICL of the model as described by Sottosanti and Risso (2021+).

ICL <- function(x){
    K <- nrow(x$mu)
    R <- ncol(x$mu)
    n <- nrow(x$x)
    p <- ncol(x$x)
    max(x$logL) - n*log(K) - p*log(R) - .5*(4*K*R+R)*log(n * p)
}
