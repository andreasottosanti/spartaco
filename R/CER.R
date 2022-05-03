#' Classification Error Rate (CER)
#'
#' This function returns the classification error rate (CER, *Witten and Tibshirani, 2010*) as a measure of the discrepancy between a reference classification and an estimated classification.
#'
#' @export
#'
#'
#'@param reference the vector containing the reference classification.
#'@param estimate the vector containing the estimated classification.
#'
#'@return The function returns the CER. The closer is the value to 0,
#'the better is the estimated classification with respect to the reference one.

#' @references
#' Witten, D. M. and Tibshirani, R. (2010). A Framework for Feature Selection in Clustering. *Journal of the American Statistical Association* **105** 713-726.

#' @examples
#' a <- c(1, 0, 0, 1, 1)
#' b <- c(1, 1, 0, 1, 1)
#' CER(a, b)

CER <- function(reference, estimate){
    if(length(reference) != length(estimate)) stop("The two objects have different length")
    n <- length(reference)
    value <- 0
    for(i in 1:(n-1)){
        mP <- mQ <- numeric(n-i)
        for(j in (i+1):(n)){
            if(reference[i] == reference[j]) mP[j-i] <- 1
            if(estimate[i] == estimate[j]) mQ[j-i] <- 1
        }
        value <- value + sum(abs(mP - mQ))
    }
    return(value/(n*(n-1)/2))
}
