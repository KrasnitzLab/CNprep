#' @title Idenfication oflower threshold segments to retain
#' 
#' @description This function idenfies the segments for the lower threshold 
#' that also contain at least one segment that respect the upper threshold.
#' 
#' @param vstart a \code{vector} of \code{numeric} representing the 
#' left boundary positions of segments that respect the upper threshold limit.
#' 
#' @param vend a \code{vector} of \code{numeric} representing the 
#' right boundary positions of segments that respect the upper threshold limit.
#' 
#' @param wstart a \code{vector} of \code{numeric} representing the 
#' left boundary positions of segments that respect the lower threshold limit.
#' 
#' @param wend a \code{vector} of \code{numeric} representing the 
#' right boundary positions of segments that respect the lower threshold limit.
#' 
#' @return a \code{matrix} with a number of rows corresponding to the length of
#' \code{wstart} and with 2 columns:
#' \itemize{
#' \item {\code{startafterstart}}{ a \code{numeric} identifying the upper 
#' segment that starts after the lower segment associated to the row.}
#' \item {\code{endbeforeend}}{ a \code{numeric} identifying the upper segment 
#' that ends before the lower segment associated to the row.}
#' }
#'
#' @details TODO
#' 
#' @examples
#'
#' ## Vectors of left and right boundary positions for segments that respect
#' ## the upper threshold limit
#' upperStart <- c(75794987, 87695620, 88864215, 111800683)
#' upperEnd   <- c(75809906, 87762703, 95041220, 111898394)
#' 
#' ## Vectors of left and right boundary positions for segments that respect
#' ## the lower threshold limit
#' lowerStart <- c(75794987, 87647882, 88805625, 111799423, 153624187, 
#'                 184116712)
#' lowerEnd   <- c(75809906, 87763963, 95062284, 112097412, 153624592, 
#'                 184150457)
#' 
#' ## TODO
#' CNprep:::containment.indicator(vstart=upperStart, vend=upperEnd, 
#'                                     wstart=lowerStart, wend=lowerEnd)
#' 
#' @author Alexander Krasnitz, Guoli Sun
#' @keywords internal
containment.indicator <- function(vstart, vend, wstart, wend) 
{
    lw <- length(wstart)
    lv <- length(vstart)
    z <- cbind(c(vend, wend), c(seq_len(lv), rep(0, lw)), 
                                    c(rep(0, lv), seq_len(lw)))
    z <- z[order(z[, 1]), ]
    endbeforeend <- cummax(z[, 2])[order(z[, 3])][sort(z[, 3]) != 0]
    z <- cbind(c(wstart, vstart), c(rep((lv + 1), lw), seq_len(lv)), 
                                        c(seq_len(lw), rep(0, lv)))
    z <- z[order(z[, 1]), ]
    startafterstart <- rev(cummin(rev(z[, 2])))[order(z[, 3])][sort(z[, 3]) 
                                                                        != 0]
    return(cbind(startafterstart, endbeforeend))
}
