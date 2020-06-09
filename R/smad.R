#' @title Calculate the median absolute deviation of the values from
#' a contiguous subsection of specified vector.
#' 
#' @description Compute the median absolute deviation for a specified vector of
#' values.  Only a contiguous subsection of the vector is used for the 
#' calculation, as the first and last position are set by user.
#' 
#' @param pos a \code{vector} of 2 \code{integer} that represent the first and 
#' last positions of \code{vector} \code{v} to used for the calculation.
#' 
#' @param v a \code{vector} of \code{double} containing the values used for
#' the calculation. However, only a subsection of the \code{vector}, as set 
#' by \code{pos}, is used.
#' 
#' @param w \code{vector} of \code{double} containing the values used for
#' the weight. However, only a subsection of the \code{vector}, as set 
#' by \code{pos}, is used. When \code{NULL}, the weight is not used in
#' the calculation.
#' Default: \code{NULL}.
#' 
#' @param cN a \code{double} a scale factor for the weighted mad
#' Default: \code{1.4826}.
#' 
#' @return a \code{double} which is the median of the values.
#' 
#' @examples
#' 
#' ## A vector with the first and last positions to subset the value vector
#' position <- c(1, 8)
#' 
#' ## A value vector used to do the calculation
#' values <- c(0.172073840, 0.219913919, 0.154459489, 0.040994620, -0.182843732,
#'     0.093052725, 0.170908930, 0.086624752, -0.003855011, -0.195791649)
#'     
#' ## Calculate the median absolute deviation of the values from the 
#' ## subsetted vector
#' CNprep:::smad(pos=position, v=values)
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats mad
#' @keywords internal
smad <- function(pos, v, w=NULL, cN=1.4826) {
    res <- NULL
    # Modified for weight
    if (length(w) > 0) {
        res <- weighted.median(abs(v[pos[1]:pos[2]] -
                                       weighted.median(v[pos[1]:pos[2]], w[pos[1]:pos[2]])),
                               w[pos[1]:pos[2]]) * cN
    } else {
        res <- mad(v[pos[1]:pos[2]], na.rm = TRUE)
    }
    return(res)
}
