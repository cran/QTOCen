#' @title Function to generate simulation data from a sequentially randomized
#'  experiment designed in \insertCite{jiang2017estimation}{QTOCen}
#'
#' @details This generative model is proposed in \insertCite{jiang2017estimation}{QTOCen}, Section 5, the second example. 
#' It uniformly defined three
#' sets of conditional distributions of the survival times given the observable
#' covariates at each stage within the same framework.
#' 
#' All three models satisfy the independent censoring assumption.
#' 
#' @return 
#' This function returns a data.frame with simulated subject trajectories.
#' \itemize{
#'  \item{\code{x0}}{ the baseline covariate, always observable at relative time point 0; }
#'  \item{\code{a0}}{ the observed first-stage treatment level at relative time point 0;}
#'  \item{\code{x1}}{ an updated covariate observable to the relative time point 
#'                    \code{s_Diff_Time}, when the 
#'                    a second stage treatment is scheduled}
#'  \item{\code{a1}}{ the observed second-stage treatment level at relative time point \code{s_Diff_Time}.}
#' }
#' 
#' 
#' @param n sample size
#' @param case string. One of {"a", "b", "c"}, corresponding to three models.
#' @param s_Diff_Time Numeric. Default is 1. This is the length of time between two stages of treatment
#' @param C_max Numeric. Default is 5. This the upper bound of the uniform distribution of the 
#' censoring time variable. Changing this value shifts the overall censoring rate easily.
#' @param Censored Boolean. Default is TRUE. Whether the data has censoring or not. If TRUE, all survival time
#' would not be censored at all in the returned data.
#' @param fix_x0_value Numeric. Default is Null. If supplied, it will generate simulated 
#' data with a fixed value, \code{fix_x0_value}, of the univariate baseline covarate. 
#' 
#' 
#' @export
#' @examples
#' dataA <- simJLSDdata(500,case="a")
#' dataB <- simJLSDdata(500,case="b")
#' dataC <- simJLSDdata(500,case="c")
#' 
#' @references 
#' \insertRef{jiang2017estimation}{QTOCen}
#'
simJLSDdata <- function(n, case="a", s_Diff_Time=1, C_max = 5, Censored = TRUE,
                        fix_x0_value=NULL){
  a0 <- rbinom(n = n, 1, prob = 0.5)
  if(is.null(fix_x0_value)){ x0 <- runif(n = n, min = 0, max=4)} else {
    x0 <- fix_x0_value}
  x1 <- 0.5*x0 - 0.4*(a0 - 0.5) + runif(n = n, min = 0, max = 2)
  if(Censored) {C <- runif(n = n, min = 0, max=C_max)} else {
    C <- rep(1000, n)}
  a1 <- rbinom(n = n, 1, prob = 0.5)

  if(case == "a"){
    T1 <- rexp(n = n, rate = 0.5*exp(1.75*(a0-0.5)*(x0-2)))
    T2 <- rexp(n = n, rate = 0.3*exp(2.5*(a1-0.4)*(x1-2) - a0*(x0-2))  )

  } else if(case == "b"){
    T1 <- rexp(n = n, rate = 0.1*exp(2 *(a0-0.5)*(x0-2)))
    T2 <- rexp(n = n,      rate = 0.2*exp(3*(a1 -0.4)*(x1-2) - 3*(a0-0.5)*(x0-2))  )

  } else if(case == "c"){
    # T1 <- rexp(n = n, rate = 0.2*exp(1.5 *(a0-0.3)*(x0-3)))
    # T2 <- rexp(n = n, rate = 0.3*exp(2*(a1-0.5)*(x1-2) + 0.5*(a0-0.7)*(x0-1))  )
    T1 <- rexp(n = n, rate = 0.3*exp(3 *(a0-0.3)*(x0-3))  )
    T2 <- rexp(n = n, rate = 0.3*exp(2*(a1-0.5)*(x1-2)   - 0.5*(a0-0.3)*(x0-3))  )

  }else if(case == "b2"){
    T1 <- rexp(n = n, rate = 0.2*exp(2 *(a0-0.5)*(-x0 +2)))
    T2 <- rexp(n = n, rate = 0.2*exp(1.5 *(a1-0.5)*(x1-2)+ 0.3*x1 + 0.3*x0) )

  }

  # Z: eligibility for second stage randomization
  Z <- (pmin(C, T1) > s_Diff_Time)
  T_surv <- Z*(s_Diff_Time + T2) + (1-Z)*T1
  delta <- as.numeric(T_surv<=C)
  censor_y <- pmin(T_surv, C)
  # plot(survfit(Surv(censor_y, 1 - delta)~1))
  a1_obs <- a1; a1_obs[which(!Z)] <- NA
  x1_obs <- x1; x1_obs[which(!Z)] <- NA

  mean(delta)
  data.frame(a0 =a0, x0=x0, a1 = a1_obs, x1=x1_obs, censor_y, delta, Z,
             T1, T2,
             C, T_surv)
}


