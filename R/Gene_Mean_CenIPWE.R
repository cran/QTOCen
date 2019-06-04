#' @title A low-level function for the generic optimization step in estimating Mean-optimal
#' treatment regime for censored data
#' 
#' @description This function supports the \code{IPWE_mean_IndCen} function.
#' It does the genetic algorithm based method with inverse probability weighting for censored data. 
#' In the future, if more complicated applications/scenarios is sought after for mean optimality, 
#' users may create their own wrapper function
#' based on \code{Gene_Mean_CenIPWE}.
#' 
#' @param data_aug a data.frame of the observed data after preprocessing. It should include be
#' augmented with two new columns: \code{ph} for the enstimated propensity scores and 
#' \code{ghat} for the estimated conditional survival probabilities.
#' @param ph propensity score estimates. For example, if the treatment is denoted by \code{A},
#' then \code{ph} should be P(A=1|X)
#' @param p_level printing level
#' @param cluster default is FALSE. This can also be an object of the 'cluster' class 
#' returned by one of the makeCluster commands in the parallel package or
#'  a vector of machine names so rgenoud::genoud can setup the cluster automatically.
#' @param Domains default is NULL. Otherwise, the object should be a \code{nvars *2} 
#' matrix used as the space of parameters, which will be supplied to \code{rgenoud::genoud}.
#' @param regimeClass a formula indicating the form of treatment regimes
#' @inheritParams IPWE_Qopt_IndCen
#' @importFrom rgenoud genoud
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
#' 
#' # estimate the mean-optimal treatment regime
#' meanopt_fit <- Gene_Mean_CenIPWE(data=data_aug, ph = data_aug$ph, p_level=1, regimeClass=a~x1*x2) 
#' 
#' @export
Gene_Mean_CenIPWE <- function(data_aug, ph, p_level,
                               regimeClass,  Domains=NULL, cluster = FALSE, s.tol = 1e-04, it.num = 8,
                               pop.size = 3000) {

  regimeClass <- as.formula(regimeClass)
  txname <- as.character(regimeClass[[2]])
  DsgnMtx <- model.matrix(regimeClass, data_aug)

  nvars <- ncol(DsgnMtx)
  if (is.null(Domains)) {
    Domains <- cbind(rep(-3000, nvars), rep(3000, nvars))
  } else  {
    if (!all(dim(Domains) == c(nvars, 2)))
      stop("The dimension of customized Domains should be:\n (number of parameters-1, 2).")
  }
  
  est <- genoud(
    fn = est_mean_ipwe,
    nvars = nvars,
    x = DsgnMtx,
    censor_y = data_aug$censor_y,
    delta = data_aug$delta,
    a = data_aug$a,
    ph = ph,
    ghat = data_aug$ghat,
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
    Domains = Domains,
    starting.values = c(0.2, rep(0.5, nvars - 1)),
    hard.generation.limit = F,
    solution.tolerance = s.tol,
    optim.method = "Nelder-Mead",
    cluster = cluster
  )

  ######### estimated coefficient ####################
  horowitz_norm <- function(x) {
    if(x[2]!=0){ return(x/abs(x[2])) }else{
      stop("The estimated parameter for the first non-intercept covariate is zero.
           To satisfy the condition of horotiwz(1992) normalization method, please
           either drop this variable or change the order of variables in 'regimeClass'")
    }
  }
  
  coefficients <- est$par
  if (prod(coefficients == rep(0, nvars)) == 1)
    coefficients <- rep(0, nvars) else coefficients <- horowitz_norm(coefficients)
  hatQ <- est$value

  output <- list(coefficients = coefficients, hatQ = hatQ)
  return(output)
}
