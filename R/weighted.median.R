#' @title Calculate the lower weighted median of a specified vector.
#' 
#' @description Compute the lower weighted median of a specified vector.
#' 
#' @param v a \code{vector} of \code{double} containing the values used for
#' the calculation.
#' 
#' @param weights a \code{vector} of \code{double} containing the weight for 
#' each value of the \code{v} vector.
#' 
#' @return a \code{double} which is the weighted median of the vector.
#' 
#' @examples
#' 
#' ## Vector of values
#' values <- c(0.06363411, 0.04342896, 0.07207384, 0.07319237, 0.07273546,
#'     0.01463932, 0.02136043, 0.01967027)
#' 
#' ## Vector of weigth associated to each value    
#' weight <- c(1, 1, 1, 1, 0, 0, 0, 0)
#' 
#' ## Calculation of the weighted median
#' CNprep:::weighted.median(v = values, weights = weight)
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @keywords internal
weighted.median <- function(v, weights)
{
    weights <- weights[order(v)]
    v <- sort(v)
    sw <- sum(weights)
    return(v[which.min(abs(cumsum(weights) - 0.5 * sw))])
}
