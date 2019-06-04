#' @title Estimate the marginal quantile response of a linear static treatment regime
#' 
#' @description Assume we have binary treatment options for each subject in the target population.
#' This function evaluates a given treatment regime by the estimated 
#' marginal mean response. 
#' We assume the space of treatment regimes are linear
#' decision functions indexed by parametric coefficients.
#' 
#' This R function is an empirical \emph{value function} in the 
#' literature of optimal treatment regime estimation. Since the goal here
#' is to maximize population's \strong{marginal quantile}, this function, which estimates the perforamce
#' of a set of parameters in terms of \strong{marginal quantile},
#' is the objective function in a nonparametric policy-search method.
#' 
#' The user facing application which utilizes this function is \code{\link{IPWE_Qopt_IndCen}}.
#' 
#' @inheritParams est_mean_ipwe
#' @inheritParams IPWE_Qopt_IndCen
#' 
#' @param beta Numerical vector. Exclude the coefficient for the first nontrivial covariate. So if
#' there are \code{k} covariates, the length of \code{beta} should equal \code{k+1-1=k} because
#' the intercept needs one coefficient as well.
#' 
#' @param sign_beta1 logical. FALSE if the coefficient for the first continuous variable 
#' is fixed to be negative one; TRUE if positive one.
#' 
#' @param epsi  the product of (1) the probability of being assigned the 
#' observed treatment level through the original treatment assignment mechanism 
#' and (2) the conditional survival probability of the censoring variable at \code{censor_y}.
#' 
#' @param Penalty.level the level that determines which objective function to use. 
#' \code{Penalty.level = 0} indicates no regularization;
#' \code{Penalty.level = 1} indicates the value function estimation minus the means absolute average coefficient 
#' is the output, which is useful trick to achieve uniqueness of estimated optimal TR
#' when resolution of input response is low.
#' 
#' @importFrom quantreg rq
#' @export
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
#' # here the value for argument epsi uses 0.5 vector for brevity in notation.
#' quant_hat <- est_quant_ipwe(beta=c(-1,2), sign_beta1=TRUE, x=cbind(1, data$x1, data$x2), 
#'                             censor_y = data$censor_y, delta = data$delta, tau=0.5,
#'                             epsi = rep(0.5,n), a = data$a)  
#'
est_quant_ipwe <- function(beta,
                           sign_beta1,
                           x,
                           censor_y,
                           delta,
                           epsi,
                           a,
                           tau,
                           check_complete = TRUE,
                           Penalty.level = 0) {
  beta_full <-if (sign_beta1) append(beta, 1, after = 1) else append(beta, -1, after = 1)
  g <- as.numeric(x %*% beta_full > 0)
  nonMissing_a.d <- (a == g)
  if (check_complete) {
    if (mean(nonMissing_a.d) < 0.05) {
      message("Warning: less than 5% is consistent with observed data")
      message("A sufficiently small value is returned")
      return(-1000 * max(abs(censor_y)))
    }
  }

  R <- (a * g + (1 - a) * (1 - g)) * delta
  
  epsi[epsi < 0.02] <- 0.02
  wts <- R * (1/epsi)

  fit <- rq(censor_y ~ 1, tau, weights = wts)
  if(Penalty.level==0){
    return(as.numeric(fit$coef[1]))
  } else if(Penalty.level==1) {
    return(as.numeric(fit$coef[1]) -  mean(abs(g)))
  } else {message("wrong level of penalty")}
}
