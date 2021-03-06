#' @title Estimate the marginal mean response of a linear static treatment regime
#'
#' @description Assume we have binary treatment options for each subject in the target population.
#' This function evaluates a given treatment regime by the estimated 
#' marginal mean response. 
#' We assume the space of treatment regimes are linear
#' decision functions indexed by parametric coefficients.
#' 
#' This R function is an empirical \emph{value function} in the 
#' literature of optimal treatment regime estimation. Since the goal here
#' is to maximize population's \strong{marginal mean} response, this function, which estimates 
#' the performance of a set of parameters in terms of the \strong{marginal mean},
#' is the objective function in a nonparametric policy-search method.
#' 
#' The user facing application which utilizes this function is \code{\link{IPWE_mean_IndCen}}.
#' 
#' @param beta Numeric vector. A set of parameter that indexes the regime.
#' @param x Numeric Matrix. The baseline covariates from all observed data.
#' @param censor_y Numeric vector. The censored survival times from all observed data, i.e. \code{censor_y = min(Y, C)}
#' @param delta Numeric vector. The censoring indicators from all observed data. We use 1 for uncensored, 0 for censored.
#' @param ph Numeric vector. The estimated propensity score of being assigned treatment \code{A=1} 
#' by the original data generation mechanism for all observed data.
#' @param a Numeric vector. The vector of observed treatment level for all observed data. Treatment levels
#' should be coded as 0/1.
#' 
#' @param ghat Numeric vector. The conditional/unconditional probabilities of 
#' event that the censoring variable is larger than the observed survival time given covariates 
#' for each observation.
#' a.k.a  \eqn{F(T > y_0 \mid x_0).}
#' This can be calculated by function \code{\link{LocalKM}}. 
#' Estimation of conditional cumulative function value at \eqn{y_0} is
#' implemented in \code{\link{tauhat_func}}. 
#' 
#' @param check_complete logical. Since this value estimation method is purely
#' nonparametric, we need at least one unit in collected data such that the observed
#' treatment assignment is the same what the regime parameter suggests. If \code{check_complete}
#' is \code{TRUE}. It will check if any observation satisfies this criterion. 
#' When none observation satisfies, a message is printed to console to raise users
#' awareness that the input regime parameter \code{beta} does not agree with any observed treatment level assignment.
#' Then a sufficiently small number is returned from this function, to keep
#' the genetic algorithm running smoothly.
#' 
#' @export
#' 
#' @examples 
#' GenerateData <- function(n)
#' {
#'   x1 <- runif(n, min=-0.5,max=0.5)
#'   x2 <- runif(n, min=-0.5,max=0.5)
#'   error <- rnorm(n, sd= 1)
#'   ph <- rep(0.5,n)
#'   a <- rbinom(n = n, size = 1, prob=ph)
#'   c <- 1.5 +  + runif(n = n, min=0, max=2)
#'   cmplt_y <-  pmin(2+x1+x2 +  a*(1 - x1 - x2) +  (0.2 + a*(1+x1+x2)) * error, 4.4)
#'   censor_y <- pmin(cmplt_y, c)
#'   delta <- as.numeric(c > cmplt_y)
#'   return(data.frame(x1=x1,x2=x2,a=a, censor_y = censor_y, delta=delta))
#' }
#' n <- 100
#' data <- GenerateData(n)
#' 
#' # here the value for argument ghat uses 0.5 vector for brevity.
#' mean_hat <- est_mean_ipwe(c(-1,0,2), x=cbind(1, data$x1, data$x2), 
#'                           censor_y = data$censor_y, delta = data$delta, ph = rep(0.5,n), 
#'                           a = data$a, ghat = rep(0.5,n))  
#' 
#'
est_mean_ipwe <- function(beta, x, censor_y, delta, ph,
                           a,  ghat, check_complete=TRUE) {

  g <- as.numeric(x %*% beta > 0)
  nonMissing_a.d <- (a==g)
  if(check_complete){
    if (!any(nonMissing_a.d)) {
      message("Since one treatment regime is counterfactual
              for all observations, \n we cannot estimate
              the margianl mean of potential outcome")
      message("A sufficiently small number is returned")
      return(-1000 * max(abs(censor_y)))
    }
  }

  R <- (a * g + (1 - a) * (1 - g)) * delta
  epsi <- (ph * g + (1 - ph) * (1 - g)) * ghat

  epsi[epsi < 0.02] <- 0.02
  wts <- R * (1/epsi)

  fit <- weighted.mean(censor_y, w=wts)
  return(fit)
}
