#' @title Kernel-based Local Kaplan-Meier Estimator 
#' 
#' @description This is the local KM estimator customized for this library to run
#' in batch mode.
#' It returns the estimated conditional survival probabilities given a user specified
#' set of covariate names that the survival time depends on,
#' a.k.a  \eqn{F(T > y_0 \mid x_0).}
#' 
#' More specifically, for uncensored data points, we return \code{ (1 - \link{tauhat_func}()) }.
#' If the observed data point is censored, then this function returns value -1
#' as a flag meaning we cannot .
#' 
#' 
#' @param D a data.frame with column \code{censor_y}, column \code{delta}, and additional covaraites.
#' @param NamesCov the vector of column names in data.frame \code{D} such that the
#' survival time depends on. 
#' @param bw the bandwidth parameter
#' 
#' @examples 
#' GenerateData <- function(n)
#' {
#'   x1 <- runif(n, min=-0.5,max=0.5)
#'   x2 <- runif(n, min=-0.5,max=0.5)
#'   error <- rnorm(n, sd= 1)
#'   ph <- exp(-0.5+1*(x1+x2))/(1+exp(-0.5 + 1*(x1+x2)))
#'   a <- rbinom(n = n, size = 1, prob=ph)
#'   c <- 1.5 +  + runif(n = n, min=0, max=2)
#'   cmplt_y <-  pmin(2+x1+x2 +  a*(1 - x1 - x2) +  (0.2 + a*(1+x1+x2)) * error, 4.4)
#'   censor_y <- pmin(cmplt_y, c)
#'   delta <- as.numeric(c > cmplt_y)
#'   return(data.frame(x1=x1,x2=x2,a=a, censor_y = censor_y, delta=delta))
#' }
#' n <- 20
#' D <- GenerateData(n)
#' mean_hat <- LocalKM(D, 5, c("x1","x2"))
#' 
#' 
#' @export
#' @return A vector of estimated conditional survival probability evaluated at the 
#' observed actual survival time on the same individual
#' 

LocalKM <- function(D, bw, NamesCov) {
    out <- NULL

    for (i in 1:nrow(D)) {
        x0 <- as.numeric(D[i, NamesCov])
        y0 <- D[i, "censor_y"]
        delta0 <- D[i, "delta"]
        
        if(delta0==0) {
          out <- c(out, -1)
        } else {
          # the input delta variable indicates whether the true
          # censoring time is observable
          tau.star <- tauhat_func(y0 = y0, x0 = x0, z = D$censor_y,
                                  x = as.matrix(D[, NamesCov, drop = FALSE]),
                                  delta = 1 - D$delta,
                                  bw = bw)
          out <- c(out, 1 - tau.star)  # output the estimated survival probability
        }
    }
    return(out)
}
