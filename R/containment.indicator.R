#' @title TODO
#' 
#' @description TODO
#' 
#' @param vstart TODO 
#' 
#' @param vend TODO
#' 
#' @param wstart TODO
#' 
#' @param wend TODO
#' 
#' @return TODO
#'
#' @examples
#'
#' # TODO
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
