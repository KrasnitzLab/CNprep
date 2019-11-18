#' @title Calculate the median of the values from
#' a contiguous subsection of specified vector.
#' 
#' @description Calculate the median of the values from a specified vector. 
#' Only a contiguous subsection of the vector is used for the calculation, as 
#' the first and last position are set by user.
#' 
#' @param pos a \code{vector} of 2 \code{integer} that represent the first and 
#' last positions of \code{vector} \code{v} to used for the calculation.
#' 
#' @param v a \code{vector} of \code{double} containing the values used for
#' the calculation. However, only a subsection of the \code{vector}, as set 
#' by \code{pos}, is used.
#' 
#' @return  a \code{double} which is the median of the values.
#' 
#' @examples
#' 
#' ## A vector with the first and last positions to subset the value vector
#' position <- c(1, 10)
#' 
#' ## A value vector used to do the calculation
#' values <- c(0.072073840, 0.119913919, 0.154459489, 0.040994620, -0.082843732,
#'     0.093052725, 0.170908930, 0.086624752, -0.003855011, -0.195791649,
#'     0.012126180, 0.043428961, 0.028435453, 0.075708220, 0.020358061)
#' 
#' ## Calculate the median of the values from the subsetted vector
#' CNprep:::smedian(pos=position, v=values)
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats median
#' @keywords internal
smedian <- function(pos, v) {
    median(v[pos[1]:pos[2]], na.rm=TRUE)
}
