#' SpaRTaCo ICL
#'
#' This function returns the integrated completed log-likelihood (ICL, *Biernacki et al., 2000*) of a SpaRTaCo model.
#'
# #' @import SpatialExperiment
#' @export
#'
#' @param x an object of class `spartaco`.
#'
#' @return The ICL of the model as described by *Sottosanti and Risso (2021+)*.
#'
#' @references
#' Biernacki, C., Celeux, G. and Govaert, G. (2000). Assessing a mixture model for clustering with the
#' integrated completed likelihood. *IEEE transactions on pattern analysis and machine intelligence* **22** 719–725.
#'
#'  Sottosanti, A. and Risso, D. (2021+) Co-clustering of Spatially Resolved Transcriptomic Data [(preprint)](https://arxiv.org/abs/2110.04872)

ICL <- function(x){
    K <- nrow(x$mu)
    R <- ncol(x$mu)
    n <- nrow(x$x)
    p <- ncol(x$x)
    if(is.null(x$logL)) max(x$max.logL) - n*log(K) - p*log(R) - .5*(4*K*R+R)*log(n * p) else
        max(x$logL) - n*log(K) - p*log(R) - .5*(4*K*R+R)*log(n * p)
}
