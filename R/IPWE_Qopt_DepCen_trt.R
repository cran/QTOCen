#' @title Estimate the Quantile-opt Treatment Regime under the assumption that the censoring
#' time's distribution only depends on treatment level
#'
#' @description   Here we assume the censoring variable is independent of covariates 
#' and potential outcomes 
#' given the treatment assignment. For example, if evidence shows that patients at 
#' certain treatment level are prone to experience censoring earlier.
#' 
#' @param data raw data.frame
#' @param tau the quantile of interest
#' @inheritParams IPWE_Qopt_IndCen
#'
#' @importFrom rgenoud genoud
#' @importFrom survival Surv survfit
#' @export
#' @details  data is a dataframe that contains:
#'     a(observed treatment assignment), 
#'     censor_y, and
#'     delta
#'  
#' @examples
#' GenerateData_DepCen_trt <- function(n)
#' {
#'   x1 <- runif(n, min=-0.5,max=0.5)
#'   x2 <- runif(n, min=-0.5,max=0.5)
#'   error <- rnorm(n, sd= 1)
#'   ph <- exp(-0.5+1*(x1+x2))/(1+exp(-0.5 + 1*(x1+x2)))
#'   a <- rbinom(n = n, size = 1, prob=ph)
#'   c <- 1 + 1*a + runif(n = n, min=0, max=2)
#'    # distribution of `c' depends on treatment level `a'
#'   cmplt_y <-  pmin(2+x1+x2 +  a*(1 - x1 - x2) +  (0.2 + a*(1+x1+x2)) * error, 4.4)
#'   censor_y <- pmin(cmplt_y, c)
#'   delta <- as.numeric(c > cmplt_y)
#'   return(data.frame(x1=x1,x2=x2,a=a, censor_y = censor_y, delta=delta))
#' }
#' \dontshow{
#' data <- GenerateData_DepCen_trt(50)
#' fit2 <- IPWE_Qopt_DepCen_trt(data = data, regimeClass = a~x1+x2, moPropen = a~x1+x2,
#'                                  tau = 0.2, pop.size=300, it.num = 3)
#' }                        
#' \donttest{
#' n <- 400
#' data <- GenerateData_DepCen_trt(n)
#' fit1 <- IPWE_Qopt_DepCen_trt(data = data, regimeClass = a~x1+x2, moPropen = a~x1+x2,
#'                                  tau = 0.2)
#'                                  }
#'    

IPWE_Qopt_DepCen_trt <- function(data, regimeClass, tau,
                            moPropen = "BinaryRandom",
                            cluster = FALSE, p_level = 1, s.tol = 1e-04, it.num = 8,
                            pop.size = 6000) {

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
    yNA.idx <- which(is.na(data$y))
    data <- data[!is.na(data$censor_y), ]
    message(paste("(", numNAy, "observations are removed since outcome is missing)"))
  }


  n <- nrow(data)
  regimeClass <- as.formula(regimeClass)
  txname <- as.character(regimeClass[[2]])
  if (length(txname) != 1)
    stop("Error: only one term should be used to represent the treatment")
  if (!all(unique(data[, txname]) %in% c(0,1)))
    stop("Error: the levels of treatment should be binary numeric: 0 and 1.")

  # Propensity score
  if (moPropen == "BinaryRandom") {
    data$ph <- rep(mean(data$a), times=n)
  } else {
    moPropen <- as.formula(moPropen)
    logistic.model.tx <- glm(moPropen, data = data, family = binomial)
    data$ph <- as.vector(logistic.model.tx$fit)
  }


  D_0 <- data[data[, txname] == 0, ]
  D_1 <- data[data[, txname] == 1, ]
  # KM estimates of distribution of the censoring variable for the
  # two treatment groups.
  # Notice that the censoring variable is missing if delta ==1.
  # survfit(Surv(censor_y, event = deltaC)~1, data=data)
  # survest <- stepfun(survfit_all$time, c(1, survfit_all$surv))
  # data$ghat <- survest(data$censor_y)
  survfit_0 <- survfit(Surv(censor_y, 1 - delta)~1, data=D_0)
  survfit_1 <- survfit(Surv(censor_y, 1 - delta)~1, data=D_1)
  survest_0 <- stepfun(survfit_0$time, c(1, survfit_0$surv))
  survest_1 <- stepfun(survfit_1$time, c(1, survfit_1$surv))

  D_0$ghat <- survest_0(D_0$censor_y)
  D_1$ghat <- survest_1(D_1$censor_y)

  # add either the true censoring probability G, or an
  # estimation
  pD <- rbind(D_0, D_1)
  pD$epsi <- (pD$ph * pD[, txname] + (1 - pD$ph) * (1 - pD[, txname])) * pD$ghat
  
  fit_tau <- Gene_Quantile_CenIPWE(data_aug = pD,  
                                   tau = tau, regimeClass = regimeClass, cluster = cluster,
                                   pop.size = pop.size, p_level = p_level, s.tol = s.tol,
                                   it.num = it.num)
  fit <- NULL
  
  # Horowitz reparametrization
  coefficient <- fit_tau$coefficient / abs(fit_tau$coefficient[2])
  names(coefficient) <- colnames(model.matrix(regimeClass, data))
  fit$coefficients <- coefficient
  fit$hatQ <- fit_tau$hatQ
  fit$moPropen <- moPropen
  fit$regimeClass <- regimeClass
  fit$tau <- tau
  fit$data_aug <- data

  class(fit) <- c("Qopt", "Censored", "Treatment_dependent Censoring")
  return(fit)
}
