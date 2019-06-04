#' @title  Estimate Quantile-optimal Treatment Regime for covariates-dependent random censoring data
#' 
#' @description This function estimates the Quantile-optimal Treatment Regime 
#' for a given quantile level of interest
#' under the
#' assumption that the distribution of censoring time is independent of the set of potential
#' survival times given a set of baseline covariates and treatment actually received.
#' 
#' More specifically, we do stratification by treatment first and then
#'  used kernel smoothing to estimate local survival function of censoring time
#'  for each treatment group.
#' 
#' 
#' @param data raw data.frame
#' @param tau the quantile of interest
#' @param regimeClass the class of treatment regimes. e.g., 'txname ~ x1+x2'.
#' @param moPropen an optional string for the working model of treatment assignment
#' @param DepCens an optional vector of baseline variable names that the censoring variable 
#' depends on. Note that
#' the treatment variable is always treated as dependent with the censoring time. 
#' If unspecified (\code{DepCens=NULL}), then all variables on the right side of 
#' \code{regimeClass} are used for \code{DepCens}
#' @param UseTrueG logical. Whether the true survival probability of each patient is provided. 
#' @param bw the bandwidth of local KM model (e.g. see Wang-wang 2008)
#' @param trueG_value default is NULL. 
#' IF \code{UseTrueG=FALSE}, \code{trueG_value} should be \code{NULL}.
#' @param cluster default is FALSE. This can also be an object of the 'cluster' class 
#' returned by one of the makeCluster commands in the parallel package or
#'  a vector of machine names so rgenoud::genoud can setup the cluster automatically.
#' @param p_level print level
#' @param s.tol tolerance level
#' @param it.num the maximum iteration number
#' @param pop.size the initial population size
#' @param Domains default is NULL.
#'
#' @importFrom rgenoud genoud
#' @importFrom methods is
#' @import stats
#' @export
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

#' 
#' \dontshow{
#' data <- GenerateData(50)
#' fit2 <- IPWE_Qopt_DepCen_general(data = data, regimeClass = a~x1+x2, moPropen = a~x1+x2,
#'                                  tau = 0.2, bw = 20/50, 
#'                                  pop.size=300, it.num = 3)
#' }                        
#'           
#' \donttest{
#' n <- 400
#' data <- GenerateData(n)
#' fit1 <- IPWE_Qopt_DepCen_general(data = data, regimeClass = a~x1+x2, moPropen = a~x1+x2,
#'                                  tau = 0.2, bw = 20/n, 
#'                                  pop.size=3000, it.num = 3)
#'                                  }



IPWE_Qopt_DepCen_general <- function(data,
                            regimeClass,
                            tau,
                            Domains = NULL,
                            bw,
                            moPropen = "BinaryRandom",
                            DepCens = NULL,
                            UseTrueG = FALSE,
                            trueG_value = NULL,
                            cluster = FALSE,
                            p_level = 1,
                            s.tol = 1e-05,
                            it.num = 8,
                            pop.size = 5000)
{
  if (!(exists("data") && is.data.frame(data)))
        stop("Error: data has to be a data frame")

  if (!("censor_y" %in% names(data)))
    stop("The response variable 'censor_y' is not found.")

  if (!("delta" %in% names(data)))
    stop("The censoring indicator variable 'delta' is not found.")

  if (UseTrueG) {
    if (!exists("trueG_value"))
      stop("When 'UseTrueG=TRUE', it is required to provide the value of trueG_value for each observation.")
    # augment the user-specified true survival probabilities
    data$trueG_value <- trueG_value
  }

  if (tau > 0.99 | tau < 0.01)
    stop("tau value is required to be strictly between 0.01 and 0.99")


  numNAy <- sum(is.na(data$censor_y))
  if (numNAy>0){
    yNA.idx <- which(is.na(data$censor_y))
    data<-data[!is.na(data$censor_y),]
    message(paste("(", numNAy,
                  "observations are removed since outcome is missing)"))
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
    data$ph <- rep(mean(data[,txname]), n)
  } else {
    moPropen <- as.formula(moPropen)
    logistic.model.tx <- glm(moPropen, data = data, family = binomial(link='logit'))
    data$ph <- as.vector(logistic.model.tx$fit)
  }

  if(UseTrueG) {
    data$ghat <- data$trueG_value
  } else {
    split.data <- split(data, txVec)
    for (i in 1:2) {
      x_DepCen_scaled <-
        apply(
          split.data[[i]][, all.vars(regimeClass)[-1], drop = FALSE],
          MARGIN = 2,
          FUN = function(x)
            (x - min(x)) / (max(x) - min(x))
        )
      x_DepCen_scaled2 <- as.data.frame(
        cbind(
          x_DepCen_scaled,
          censor_y = split.data[[i]]$censor_y,
          delta = split.data[[i]]$delta
        )
      )
      split.data[[i]]$ghat <- LocalKM(x_DepCen_scaled2,
                                   bw = bw,
                                   NamesCov = if(is.null(DepCens)) all.vars(regimeClass)[-1] else DepCens)
  }
  data <- rbind(split.data[[1]], split.data[[2]])
} # ends if(UsetrueG)

  txVec_reorder <- try(data[, txname], silent = TRUE)
  data$epsi <- (data$ph * txVec_reorder + (1 - data$ph) * (1 - txVec_reorder)) * data$ghat

  fit_tau <- Gene_Quantile_CenIPWE(
    data_aug = data,
    tau = tau,
    regimeClass = regimeClass,
    Domains = Domains,
    cluster = cluster,
    pop.size = pop.size,
    p_level = p_level,
    s.tol = s.tol,
    it.num = it.num
  )

  fit <- NULL
  fit$coefficients <- fit_tau$coefficients
  names(fit$coefficients) <- colnames(model.matrix(regimeClass, data))
  fit$Match_prop <- fit_tau$Match_prop
  fit$hatQ <- fit_tau$hatQ
  fit$moPropen <- moPropen
  fit$regimeClass <- regimeClass
  fit$tau <- tau
  fit$bw <- bw
  fit$data_aug <- data

  class(fit) <- c("Qopt", "Censored", "Dependent Censoring")
  return(fit)
}
