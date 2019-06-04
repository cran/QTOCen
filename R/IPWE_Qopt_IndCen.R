#' @title Function to estimate the quantile-optimal treatment regime: 
#' the independent censoring Case
#'
#' @description This function implements the estimation method proposed in Chapter 2 of
#' \insertCite{zhou2018quantile}{QTOCen}. It estimates the quantile-optimal treatment regime 
#' for a given quantile level of interest from a single-stage clinical randomized
#' experiment or 
#' a single-stage observational study under the independent censoring assumption. In other
#' words, we estimate the parameters indexing the quantile-optimal treatment regime.
#' 
#' Our assumption of independent censoring means
#' the distribution of the censoring time is the same 
#' conditional on baseline covariates, treatment group and the two potential survival times. 
#' 
#' @inheritParams IPWE_mean_IndCen
#' @param tau 	a value between 0 and 1. This is the quantile of interest.
#' @param Domains default is NULL. Otherwise, the object should be a \code{nvars *2} 
#' matrix used as the space of parameters, which will be supplied to \code{rgenoud::genoud}.
#' @param sign_beta1 logical. Default is NULL. FALSE if the coefficient for the first continuous variable is fixed to be negative one; TRUE if positive one.
#' @param Penalty.level 0: stop if the marginal quantiles cannot be further optimized; 1: continue
#' the search among treatment regimes with with same value for the TR
#' with the smallest intended proportion of treatment.
#'
#' @import survival
#' @import grDevices 
#' @importFrom rgenoud genoud
#' @export
#' @details  The input argument \code{data} is the dataframe that contains:
#' \enumerate{    
#'    \item \code{a} observed treatment assignment
#'    \item \code{censor_y} the censored response variable
#'    \item \code{delta} the censoring indicator
#'  }
#' The naming of these three columns should be strict.
#' 
#' Note that this function currently only works for scenarios in which 
#' treatment is binary.
#' 
#' @examples      
#' GenerateData <- function(n)
#' {
#'   x1 <- runif(n, min=-0.5,max=0.5)
#'   x2 <- runif(n, min=-0.5,max=0.5)
#'   error <- rnorm(n, sd= 1)
#'   ph <- exp(-0.5+1*(x1+x2))/(1+exp(-0.5 + 1*(x1+x2)))
#'   a <- rbinom(n = n, size = 1, prob=ph)
#'   c <- 1 + 1*a + runif(n = n, min=0, max=2)
#'   cmplt_y <-  pmin(2+x1+x2 +  a*(1 - x1 - x2) +  (0.2 + a*(1+x1+x2)) * error, 4.4)
#'   censor_y <- pmin(cmplt_y, c)
#'   delta <- as.numeric(c > cmplt_y)
#'   return(data.frame(x1=x1,x2=x2,a=a, censor_y = censor_y, delta=delta))
#' }
#' n <- 400
#' \donttest{
#' data <- GenerateData(n)
#' fit1 <- IPWE_Qopt_IndCen(data = data, regimeClass = a~x1+x2, tau=0.25)
#' 
#' # We can used the returned model to visualize the Kaplan-meier
#' # estimate of survival function of the censoring time variable,
#' # justified by the independent censoring assumption.
#' library(survminer)
#' ggsurvplot(fit1$survfitCensorTime, data=fit1$data_aug, risk.table = TRUE)
#'  }
#' \dontshow{
#' sdata <- GenerateData(100)
#' fit1 <- IPWE_Qopt_IndCen(data = sdata, regimeClass = a~x1, 
#'                          tau=0.25,
#'                          pop.size=500, it.num = 2, s.tol=0.5)
#'  }
#'  
#' @references
#' \insertRef{zhou2018quantile}{QTOCen}
#' 
#' \insertRef{horowitz1992smoothed}{QTOCen}

IPWE_Qopt_IndCen <- function(data,
                                regimeClass,
                                tau,
                                moPropen = "BinaryRandom",
                                Domains = NULL,
                                cluster = FALSE,
                                p_level = 1,
                                s.tol = 1e-04,
                                it.num = 8,
                                pop.size = 6000,
                                sign_beta1 = NULL,
                                Penalty.level = 0) {
  if (!(exists("data") && is.data.frame(data)))
    stop("Error: data has to be a data frame")

  if (!("censor_y" %in% names(data)))
    stop("The response variable 'censor_y' is not found in the input data.")

  if (!("delta" %in% names(data)))
    stop("The censoring indicator variable 'delta' is not found.")

  if (tau > 0.99 | tau < 0.01)
    stop("tau value is required to be strictly between 0.01 and 0.99")

  numNAy <- sum(is.na(data$censor_y))
  if (numNAy > 0) {
    data <- data[!is.na(data$censor_y), ]
    message(paste("(", numNAy, "observations are removed since outcome is missing)"))
  }

  n <- nrow(data)
  regimeClass <- as.formula(regimeClass)
  txname <- as.character(regimeClass[[2]])
  txVec <- try(data[, txname], silent = TRUE)
  if (is(txVec, "try-error")) {
    stop("Variable '", paste0(txname, "' not found in 'data'."))
  }
  if(!all(unique(txVec) %in% c(0,1)))
    stop("The levels of treatment must be numeric, being either 0 or 1.")


  # Propensity score
  if (moPropen == "BinaryRandom") {
    data$ph <- rep(mean(txVec), n)
  } else {
    moPropen <- as.formula(moPropen)
    logistic.model.tx <- glm(moPropen, data = data, family = binomial)
    data$ph <- as.vector(logistic.model.tx$fit)
  }


  # KM estimates of distribution of the censoring variable, which should be random, i.e.
  # independent of A,X
  # Notice that the censoring variable is missing if delta ==1.
  data$deltaC <- 1 - data$delta
  survfit_all <- survival::survfit(survival::Surv(time=censor_y, event = deltaC)~1, data=data)
  survest <- stepfun(survfit_all$time, c(1, survfit_all$surv))
  
  #predict Pr(C>Y_censor), to get data$ghat.
  data$ghat <- survest(data$censor_y)

  data$epsi <- (data$ph * txVec + (1 - data$ph) * (1 - txVec)) * data$ghat


  fit <- Gene_Quantile_CenIPWE(
    data_aug = data,
    tau = tau,
    regimeClass = regimeClass,
    Domains = Domains,
    cluster = cluster,
    pop.size = pop.size,
    p_level = p_level,
    s.tol = s.tol,
    it.num = it.num,
    sign_beta1 = sign_beta1,
    Penalty.level = Penalty.level
  )

  names(fit$coefficients) <- colnames(model.matrix(regimeClass, data))
  fit$moPropen <- moPropen
  fit$regimeClass <- regimeClass
  fit$tau <- tau
  fit$data_aug <- data
  fit$sign_beta1 <- sign_beta1
  fit$Penalty.level <- Penalty.level
  fit$survfitCensorTime <- survfit_all

  class(fit) <- c("Qopt", "Censored", "Independent_Censoring")
  return(fit)
}
