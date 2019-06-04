#' @title Generate biquadratic kernel weights for a univariate variable
#' @description This is the biquadratic kernel function, that weights
#' observations by their distances to the target observation.
#' @param x0k Numeric scalar. One univariate covariate value of interest 
#'            from one observation.
#' @param Xk Numerical vector. The vector of the same covariate from observations
#' @param bw.bnk The bandwith scalar parameter.
#' 
#' 
#' @note 
#' This function is widely used for generating kernel weights for 
#' nonparametrically estimating conditional survival functions. See 
#' Section 2.3 of \insertCite{wang2009locally}{QTOCen}.
#' @examples  
#' Bnk_func(x0k=0, Xk=c(-5:5), bw.bnk=10)
#' 
#' @export
#' 
#' @return This function returns a list of kernel weights with the same length of input \code{Xk}.
#' 
#' @references
#' \insertRef{wang2009locally}{QTOCen}
Bnk_func <- function(x0k, Xk, bw.bnk) {
    xx <- (Xk - x0k)/bw.bnk
    xx[abs(xx) >= 1] <- 1
    w <- 15 * (1 - xx^2)^2/16  #biquadratic kernel
    w <- w/sum(w)
    return(w)
}


