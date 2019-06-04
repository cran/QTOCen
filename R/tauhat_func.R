#' @title Kernel-based Local Kaplan-Meier Estimator for the Conditional Probability of the Survival Time
#' @importFrom Rdpack reprompt
#' 
#' @description This function estimates the value of
#'    \deqn{F(T <= y_0 \mid x_0),}
#'  the conditional cumulative distribution function of a survival time \eqn{T} 
#'  given covaraites vector \eqn{x_0} 
#'  at value \eqn{y_0}.
#' This estimator is described in detail in \insertCite{wang2009locally}{QTOCen}.
#' 
#'      
#' @param y0 the vector of censored outcome of a single observation
#' @param x0 the vector of given covariate of a single observation 
#' @param z observed vector of response variable from observed data
#' @param x the observed matrix of covariates, the dimension is # of observations by number of covariates. 
#'          Note that the vector of ones should NOT be included in \code{x}.
#' @param delta the vector of censoring indicators
#' @param bw the scalar bandwidth parameter in kernel
#' 
#' @details 
#' For cases with multivariate covariates, we adopted a product kernel. 
#' For example, in the bivariate case we use \deqn{K(x_1, x_2) = K_1(x_1) K_2(x_2),}
#' where \eqn{K_1} and \eqn{K_2} are both biquadratickernel functions.
#'  
#' 
#' @export
#' @examples 
#' tauhat_func(y0=10, x0=c(2,3), z=c(10, 12, 11), 
#'             x=matrix(c(1,1,2,2,3,3), nrow=3, byrow=TRUE), 
#'             delta=c(1,1,0), bw=10)
#' 
#' @references
#' \insertRef{wang2009locally}{QTOCen}

tauhat_func <- function(y0, x0, z, x, delta, bw) {
    # tau0(y0, x0) = F(T<y0|x0);
    # so y0 is the C_i, and x0 is the xi in the paper,
    # z is observed vector of response
    # variable x is the observed covariate delta is the
    # censoring indicator function bw.bandwidth is the
    # bandwidth
    n <- length(z)
    if (n == 0) {
        message("no input for tauhat_func")
        return(NULL)
    }

    x <- as.matrix(x)
    K <- ncol(x)
    if (K != length(x0))
      stop("Error in tauhat.func: length of x0 not equal to number of columns in x")
    Bn.array <- matrix(0, nrow = nrow(x), ncol = ncol(x))
    for (k in 1:K) {
      Bn.array[, k] <- Bnk_func(x0k = x0[k],
                                Xk = x[, k], bw.bnk = bw)  ## kernel for the k-th dimension
    }
    Bn <- apply(Bn.array, MARGIN = 1, FUN = prod)

    if (y0 < max(z)) {
        z2 = sort(z)
        Order = order(z)  # so z[Order] = z2
        Bn2 = Bn[Order]
        delta2 = delta[Order]
        eta = which(delta2 == 1 & z2 <= y0)  # the index of those observations satisfying delta2==1 & z2<=y0
        Bn3 = Bn2[n:1]  # change the order of Bn2, make the first obs of Bn2 to be the last of Bn3
        tmp = 1 - Bn2/cumsum(Bn3)[n:1]
        out = 1 - prod(tmp[eta], na.rm = T)  # na.rm=T, as some of those tmp=NA as the denom =0
    } else { out <- 1}
    return(out)
}
