#' @title A low-level function for the generic optimization step in 
#' estimating Quanilte-optimal treatment regime for censored data
#' 
#' @description This function supports several user facing functions for Quantile-optimal
#' treatment regime estimation, namely
#' 
#'  \code{\link{IPWE_Qopt_IndCen}(), \link{IPWE_Qopt_DTR_IndCen}(), \link{IPWE_Qopt_DepCen_trt}(), and \link{IPWE_Qopt_DepCen_general}()}. 
#'  
#' It implements the genetic algorithm based policy-search method with 
#' inverse probability weighting for censored data, such that the estimator is cube root consistent
#' under the assumption that the propensity score model and the model for the 
#' survival distriution of the censoring time variable are both correct.
#'  
#' @param data_aug a data.frame of the observed data after preprocessing. It should include be
#' augmented with a new column: \code{epsi} for the composite weights.
#'
#' @inheritParams est_quant_ipwe
#' @inheritParams IPWE_Qopt_IndCen
#'
#' @import survival
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
#' # preprocessing
#' data_aug <- data
#' data_aug$ph <- rep(mean(data$a), n)
#' data_aug$deltaC <- 1 - data_aug$delta
#' library(survival)
#' survfit_all <- survfit(Surv(censor_y, event = deltaC)~1, data=data_aug)
#' survest <- stepfun(survfit_all$time, c(1, survfit_all$surv))
#' data_aug$ghat <- survest(data_aug$censor_y)
#' data_aug$epsi <- (data_aug$ph * data_aug$a + (1 - data_aug$ph) * (1 - data_aug$a)) * data_aug$ghat
#' 
#' # estimate the median-optimal treatment regime
#' \donttest{
#' quantopt_fit <- Gene_Quantile_CenIPWE(data_aug=data_aug,tau=0.5,
#'                                       p_level=1, regimeClass=a~x1+x2^2, 
#'                                       sign_beta1=FALSE)
#' }
#' \dontshow{
#' quantopt_fit <- Gene_Quantile_CenIPWE(data_aug=data_aug,tau=0.5, sign_beta1=FALSE,
#'                                       p_level=1, regimeClass=a~x1+x2,
#'                                       s.tol=0.4, it.num=1, pop.size=1000)
#' }
#' @export
#' 
Gene_Quantile_CenIPWE <- function(data_aug,
                               tau,
                               p_level,
                               regimeClass,
                               cluster = FALSE,
                               s.tol = 1e-04,
                               it.num = 8,
                               pop.size = 5000,
                               Domains = NULL,
                               sign_beta1 = NULL,
                               Penalty.level = 0)
{
  regimeClass <- as.formula(regimeClass)
  txname <- as.character(regimeClass[[2]])
  DsgnMtx <- model.matrix(regimeClass, data_aug)

  # since quantile estimation is more complicated than mean estimation, two optimizatino domains are set up
  nvars <- ncol(DsgnMtx) - 1


  if (is.null(Domains)) {
    Domains <- cbind(rep(-3000, nvars), rep(3000, nvars))
  } else  {
    if (!all(dim(Domains) == c(nvars, 2)))
      stop("The dimension of customized Domains should be:\n (number of parameters-1, 2).")
  }

  est.wrapper <- function(sign_beta1) {
    est <- genoud(
      fn = est_quant_ipwe,
      sign_beta1 = sign_beta1,
      Penalty.level = Penalty.level,
      nvars = nvars,
      Domains = Domains,
      x = DsgnMtx,
      censor_y = data_aug$censor_y,
      delta = data_aug$delta,
      a = data_aug[, txname],
      epsi = data_aug$epsi,
      tau = tau,
      print.level = p_level,
      max = TRUE,
      pop.size = pop.size,
      wait.generations = it.num,
      gradient.check = FALSE,
      BFGS = FALSE,
      P1 = 50,
      P2 = 50,
      P3 = 10,
      P4 = 50,
      P5 = 50,
      P6 = 50,
      P7 = 50,
      P8 = 50,
      P9 = 0,
      P9mix = NULL,
      starting.values = NULL,
      hard.generation.limit = F,
      solution.tolerance = s.tol,
      optim.method = "Nelder-Mead",
      cluster = cluster
    )
  }

  if (!is.null(sign_beta1)) {
      est <- est.wrapper(sign_beta1 = sign_beta1)
      if (sign_beta1) {
        est$coef <- append(est$par, 1, after = 1)
      } else {
        est$coef <- append(est$par, -1, after = 1)
      }
    } else {
      est.negative <- est.wrapper(sign_beta1 = FALSE)
      est.posi <- est.wrapper(sign_beta1 = TRUE)
      if (est.negative$value > est.posi$value) {
        est <- est.negative
        est$coef <- append(est$par, -1, after = 1)
      } else {
        est <- est.posi
        est$coef <- append(est$par, 1, after = 1)
      }
    }

    ######### Output the estimation ####################
    if (sum(abs(est$coef) == 3000) > 0)
        message("Warning: the largest estimates is achieved on the boundary of the default domain")
    hatQ <- est$value
    Match_prop <- mean(as.numeric(DsgnMtx %*% est$coef > 0) == data_aug[,txname])
    Treated_prop <- mean(as.numeric(DsgnMtx %*% est$coef > 0) )

    output <- list(coefficients = est$coef, hatQ = hatQ, Match_prop=Match_prop, Treated_prop=Treated_prop)
    return(output)
}



