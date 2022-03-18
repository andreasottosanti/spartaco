#' Classification Error Rate (CER)
#'
#' This function returns the classification error rate (CER) as a measure of the discrepancy between the reference and the estimated classifications.
#'
#' @export
#'
#'
#'@param reference the vector containing the reference classification.
#'@param estimate the vector containing the estimated classification.
#'
#'@return The function returns the CER. The closer is the value to 0,
#'the better is the estimated classification with respect to the reference one.

#CER <- function(reference, estimate){
#    mP <- mQ <- matrix(0, length(reference), length(reference))
#    for(i in 1:(length(reference)-1)){
#        for(j in (i+1):length(reference)){
#            if(reference[i] == reference[j]) mP[i,j] <- mP[j,i] <- 1
#            if(estimate[i] == estimate[j]) mQ[i,j] <- mQ[j,i] <- 1
#        }
#    }
#    return(sum(abs(mP[lower.tri(mP)]-mQ[lower.tri(mQ)]))/(length(reference)*(length(reference)-1)/2))
#}

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
