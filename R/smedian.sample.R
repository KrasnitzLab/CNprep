#' @title TODO
#' 
#' @description TODO
#' 
#' @param pos TODO
#' 
#' @param v TODO
#' 
#' @return TODO
#' 
#' @examples
#' 
#' position <- c(1, 10)
#' names(position) <- c("StartProbe", "EndProbe")
#' 
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
