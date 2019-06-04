#' @title Estimate the mean-optimal treatment regime for data with independently censored response
#' @description This function estimates the Mean-optimal Treatment Regime 
#' with censored response. 
#' The implemented function only works for scenarios in which 
#' treatment is binary and the censoring time
#' is independent of baseline covariates, treatment group and all potential survival times. 
#' 
#' @param data a data.frame, containing variables in the \code{moPropen} and 
#' \code{RegimeClass} and also the response variables, namely \code{censor_y} as the censored response, 
#' and \code{delta} as the censoring indicator.
#' 
#' @param regimeClass a formula specifying the class of treatment regimes to search,
#'        e.g. if \code{regimeClass = a~x1+x2}, and then this function will 
#'        search the class of treatment regimes
#'        of the form 
#'        \deqn{d(x) = I \left(\beta_0 +\beta_1  x_1 + \beta_2  x_2 > 0\right).
#'        }{d(x)=I(\beta_0 +\beta_1 * x1  + \beta_2 * x2 > 0).}
#'        Polynomial arguments are also supported.
#'          
#'          
#' @param moPropen  The propensity score model for the probability of receiving 
#'        treatment level 1.
#'        When \code{moPropen} equals the string "BinaryRandom",  the proportion of observations
#'        receiving treatment level 1 in the sample will be plugged in as an estimate
#'        of the propensity.
#'        Otherwise, this argument should be a formula/string, based on which this function
#'        will fit a logistic regression on the treatment level.  e.g. \code{a1~x1}.          
#'          
#'          
#' @param cluster default is FALSE, meaning do not use parallel computing for the genetic algorithm(GA).
#' 
#' @param p_level choose between 0,1,2,3 to indicate different levels of output
#'          from the genetic function. Specifically, 0 (minimal printing),
#'            1 (normal), 2 (detailed), and 3 (debug).
#' 
#' @param s.tol tolerance level for the GA algorithm. This is input for parameter \code{solution.tolerance}
#' in function \code{rgenoud::genoud}.
#' 
#' @param it.num the maximum GA iteration number
#' 
#' @param pop.size an integer with the default set to be 3000. This is roughly the 
#'                number individuals for the first generation
#'                in the genetic algorithm (\code{rgenoud::genoud}).
#'                
#' @param Domains default is NULL. Otherwise, the object should be a \code{nvars *2} 
#' matrix used as the space of parameters, which will be supplied to \code{rgenoud::genoud}. 
#' \code{nvars} is the total number of parameters.
#' 
#' 
#' 
#' 
#' @return This function returns an object with 6 objects: 
#' \itemize{
#'  \item{\code{coefficients}}{ the estimated parameter indexing the mean-optimal treatment regime. 
#'          Since we focus the space of linear treatment regimes, the estimated decision rule
#'          cannot be uniquely identified without scale normalized. In this package,
#'          we normalized by \eqn{|\beta_1| = 1}, which was proposed in Horowitz \insertCite{horowitz1992smoothed}{QTOCen}.     }
#'  \item{\code{hatQ}} { the estimated optimal marginal mean response} 
#'  \item{\code{moPropen}}{ log of the input argument of \code{moPropen}}
#'  \item{\code{regimeClass}}{ log of the input argument of \code{regimeClass}}
#'  \item{\code{data_aug}}{ Training data with additional columns used in the algorithm. Note that \code{data_aug} is used for plotting 
#'  of survival function of the censoring time}
#'  \item{\code{survfitCensorTime}}{ the estimated survival function of the censoring time}
#' }
#' 
#' 
#' 
#' @importFrom rgenoud genoud
#' @importFrom survival survfit
#' @importFrom methods is
#' @import grDevices
#' @export
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
#' n <- 400
#' 
#' D <- GenerateData(n)
#' fit1 <- IPWE_mean_IndCen(data = D, regimeClass = a~x1+x2)
#'                                  
#'    
#' @references
#' \insertRef{zhou2018quantile}{QTOCen}
#' 
#' \insertRef{horowitz1992smoothed}{QTOCen}
#'                               
IPWE_mean_IndCen <- function(data, regimeClass,
                                moPropen = "BinaryRandom",
                                Domains = NULL,
                                cluster = FALSE, p_level = 1, s.tol = 1e-04, it.num = 8,
                                pop.size = 3000) {

  # tau is between 0 and 1 regimeClass should be like
  # 'txname ~ x1+x2', etc.

  if (!(exists("data") && is.data.frame(data)))
    stop("Error: data has to be a data frame")
  if (!("censor_y" %in% names(data)))
    stop("The response variable 'censor_y' is not found in the input data.")
  if (!("delta" %in% names(data)))
    stop("The censoring indicator variable 'delta' is not found.")

  numNAy <- sum(is.na(data$censor_y))
  if (numNAy > 0) {
    yNA.idx <- which(is.na(data$y))
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
    data$ph <- rep(mean(data$a), n)
  } else {
    moPropen <- as.formula(moPropen)
    logistic.model.tx <- glm(moPropen, data = data, family = binomial)
    data$ph <- as.vector(logistic.model.tx$fit)
  }


  # KM estimates of distribution of the censoring variable, which should be random, i.e.
  # independent of A,X
  # Notice that the censoring variable is missing if delta ==1.
  data$deltaC <- 1 - data$delta
  survfit_all <- survfit(Surv(censor_y, event = deltaC)~1, data=data)
  survest <- stepfun(survfit_all$time, c(1, survfit_all$surv))
  
  #predict Pr(C>Y_censor), to get data$ghat.
  data$ghat <- survest(data$censor_y)

  # add either the true censoring probability G, or an estimation
  fit_mean <- Gene_Mean_CenIPWE(data_aug = data, ph = data$ph,
                                Domains = Domains,
                                regimeClass = regimeClass, cluster = cluster,
                                pop.size = pop.size, p_level = p_level, s.tol = s.tol,
                                it.num = it.num)
  fit <- NULL
  fit$coefficients <- fit_mean$coefficients
  names(fit$coefficients) <- colnames(model.matrix(regimeClass, data))
  fit$hatQ <- fit_mean$hatQ
  fit$moPropen <- moPropen
  fit$regimeClass <- regimeClass
  fit$data_aug <- data
  fit$survfitCensorTime <- survfit_all

  class(fit) <- c("Censored", "mean_TR")
  return(fit)
}
