#' Obtain the inverse of a matrix with some rows/columns added
#'
#' @param B a squared matrix of size `p`, representing the inverse of a matrix;
#' @param K the matrix of size `p x k` added to the former matrix;
#' @param C a `k x k` matrix added to the former matrix;
#'

AddRowColInverse <- function(B, K, C){
    u1 <- K
    u2 <- B %*% u1
    d <- solve(C - t(u1) %*% u2)
    F11.inv <- B + u2 %*% d %*% t(u2)
    add <- -B %*% u1 %*% d
    rbind(cbind(F11.inv, add),
          cbind(t(add), d))
}

#' Obtain the inverse of a matrix with some rows/columns removed
#'
#' @param B a squared matrix of size `p`, representing the inverse of a matrix;
#' @param j the indexes of the rows/columns to be removed.
#'

RemoveRowColInverse <- function(B, j){
    d <- B[j,j]
    u3 <- -B[setdiff(1:ncol(B), j), j]
    u2 <- u3 %*% solve(d)
    B[-j,-j] - u2 %*% d %*% t(u2)
}
