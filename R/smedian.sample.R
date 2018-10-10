#' @title Calculate the median of the sampled values with replacement from
#' a specified vector.
#' 
#' @description Calculate the median of the sampled values with replacement 
#' from a specified vector. Only a subsection of the vector is used for the
#' sampling, as the first and last position is set by user.
#' 
#' @param pos A \code{vector} of 2 \code{integer} that represent the first and 
#' last positions of \code{vector} \code{v} to used for the sampling step.
#' 
#' @param v A \code{vector} of \code{double} containing the values used for
#' sampling. However, only a subsection of the \code{vector}, as set 
#' by \code{pos} is used.
#' 
#' @return A \code{double} which is the median of the sampled values.
#' 
#' @examples
#' 
#' ## A vector with the first and last positions to subset the value vector
#' position <- c(1, 10)
#' 
#' ## A value vector used to do sampling with replacement
#' values <- c(0.072073840, 0.119913919, 0.154459489, 0.040994620, -0.082843732,
#'     0.093052725, 0.170908930, 0.086624752, -0.003855011, -0.195791649,
#'     0.012126180, 0.043428961, 0.028435453, 0.075708220, 0.020358061)
#' 
#' CNprep:::smedian.sample(pos=position, v=values)
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats median
#' @keywords internal
smedian.sample <- function(pos, v)
{
    ## TODO: remove NA before of after selection ??????
	w <- v[pos[1]:pos[2]][!is.na(v[pos[1]:pos[2]])]
	return(median(sample(w, length(w), replace=T), na.rm=T))
}
