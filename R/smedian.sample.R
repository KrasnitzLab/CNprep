#' @title Calculate the median of the sampled values with replacement from
#' a contiguous subsection of specified vector.
#' 
#' @description Calculate the median of the sampled values with replacement 
#' from a specified vector. Only a contiguous subsection of the vector is 
#' used for the sampling, as the first and last position are set by user.
#' 
#' @param pos a \code{vector} of 2 \code{integer} that represent the first and 
#' last positions of \code{vector} \code{v} to used for the sampling step.
#' 
#' @param v a \code{vector} of \code{double} containing the values used for
#' sampling. However, only a subsection of the \code{vector}, as set 
#' by \code{pos}, is used.
#' 
#' @return a \code{double} which is the median of the sampled values.
#' 
#' @examples
#' 
#' ## A vector with the first and last positions to subset the value vector
#' position <- c(1, 10)
#' 
#' ## A value vector used to do sampling with replacement
#' values <- c(0.072073840, 0.119913919, 0.154459489, 0.040994620, 
#'     -0.082843732, 0.093052725, 0.170908930, 0.086624752, -0.003855011, 
#'     -0.195791649, 0.012126180, 0.043428961, 0.028435453, 0.075708220, 
#'     0.020358061)
#' 
#' ## Calculate the median of the sampled values from the sampled vector
#' CNprep:::smedian.sample(pos=position, v=values)
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @importFrom stats median
#' @keywords internal
smedian.sample <- function(pos, v, w = NULL)
{
    # Modified for weight
    vP <- v[pos[1]:pos[2]][!is.na(v[pos[1]:pos[2]])]
    res <- NULL
    if (length(w) > 0) {
        wP <- w[pos[1]:pos[2]][!is.na(w[pos[1]:pos[2]])]
        swP <- sum(wP)
        sel <- sample(seq_len(length(vP)), length(vP), prob=wP/swP, 
                      replace = TRUE)
        res <- weighted.median(vP[sel], wP[sel])
    } else {
        res <- median(sample(vP, length(vP), replace = TRUE), na.rm = TRUE)
    }
    
    return(res)
}
