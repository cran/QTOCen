#' @title Function to estimate the two-stage quantile-optimal dynamic treatment 
#' regime for censored data: the independent censoring Case
#' 
#' @description  This function inplements the estimator of
#' two-stage quantile-optimal treatment regime with censored outcome 
#' by inverse probability of weighting, which is proposed in Chapter 3 of
#' \insertCite{zhou2018quantile}{QTOCen}.
#' We assume the censoring is independent of everything else, including the treatment
#' covariates, and potential outcomes.
#' 
#' Specifically, we do grid search on the sign of the coefficient for the first non-intercept variables
#'  in stage 1 and stage 2 and apply genetic algorithm on the remaining coeffients simultaneously.
#'  So if stage one has d1 covariates excluding the intercept, stage two has d2, the
#'  resulting coefficient has dimension d1+d2+2.
#'
#'
#' @inheritParams IPWE_Qopt_IndCen
#' 
#' @param regimeClass.stg1 a formula specifying the class of treatment regimes for the first stage. 
#' For details of the general formulation of a linear treatment regime
#'  see \code{regimeClass} in \code{\link{IPWE_Qopt_IndCen}}. 
#' 
#' @param regimeClass.stg2 a formula specifying the class of treatment regimes for the second stage 
#' 
#' @param moPropen1 the first stage propensity score model. Default is "BinaryRandom".
#' 
#' @param moPropen2 the second stage propensity score model. Default is "BinaryRandom".
#' 
#' @param max logical. TRUE if the goal is maximization of the quantile. FALSE is the goal is minimization of the quantile.
#' 
#' @param sign_beta1.stg1 Is sign of the coefficient for the first non-intercept 
#' variable for the first stage known? Default is NULL, meaning user does not have contraint on
#' the sign;
#' FALSE if the coefficient for the first continuous variable 
#' is fixed to be \code{-1}; TRUE if \code{1}. We can make the search space discrete because we employ
#' \eqn{|\beta_1| = 1} scale normalizaion.
#' 
#' @param sign_beta1.stg2 Default is NULL. Similar to \code{sign_beta1.stg1}.
#' 
#' @param s_Diff_Time Numeric. The fixed length of time between the first stage treatment and the
#' second stage treatment
#' 
#' @param Domains1 This is optional. If not NULL, please provide 
#' the two-column matrix for the searching range of coeffients in stage one.
#' The coefficient taking value of positive/negative one should not be included.
#' 
#' @param Domains2 This is optional. If not NULL, please provide 
#' the two-column matrix for the searching range of coeffients in stage two.
#' The coefficient taking value of positive/negative one should not be included.
#'
#' @details
#' In our setting, if a subject was censored or had experienced the event of interest
#' before \code{s_Diff_Time} units of time had elapsed after the first stage of treatment,
#' s/he would not be eligible to receive a second stage treatment.
#'
#'
#' @author Yu Zhou, \email{zhou0269@umn.edu}
#' @export
#' @importFrom rgenoud genoud
#' @import stats
#' 
#' @examples 
#' \donttest{
#' D <- simJLSDdata(400, case="a")
#' fit_2stage <-IPWE_Qopt_DTR_IndCen(data=D, tau= 0.3, regimeClass.stg1 = a0~x0,
#'                      regimeClass.stg2 = a1~x1,
#'                      sign_beta1.stg1 = FALSE,
#'                      sign_beta1.stg2 = FALSE)
#' }
#' 
#' \dontshow{
#' D <- simJLSDdata(100, case="a")
#' fit_2stage <-IPWE_Qopt_DTR_IndCen(data=D, tau= 0.3, regimeClass.stg1 = a0~x0,
#'                      regimeClass.stg2 = a1~x1,
#'                      sign_beta1.stg1 = FALSE,
#'                      sign_beta1.stg2 = FALSE,
#'                      s.tol = 0.1, it.num=2, pop.size=1000)
#' }
#' 
#' @references
#' \insertRef{zhou2018quantile}{QTOCen}



IPWE_Qopt_DTR_IndCen <- function(data,
                                    tau,
                                    regimeClass.stg1,
                                    regimeClass.stg2,
                                    s_Diff_Time = 1,
                                    moPropen1 = "BinaryRandom",
                                    moPropen2 = "BinaryRandom",
                                    sign_beta1.stg1 = NULL,
                                    sign_beta1.stg2 = NULL,
                                    Penalty.level = 0,
                                    s.tol = 1e-6,
                                    it.num = 4,
                                    max = TRUE,
                                    Domains1 = NULL,
                                    Domains2 = NULL,
                                    cluster = FALSE,
                                    p_level = 1,
                                    pop.size = 10000) {
  if (!is(data, "data.frame"))
    stop("'data' must be a data frame.")

  if (!("censor_y" %in% names(data)))
    stop("The response variable 'censor_y' is not found in the input data.")

  if (!("delta" %in% names(data)))
    stop("The censoring indicator variable 'delta' is not found.")

  if (tau > 0.99 | tau < 0.01)
    stop("tau value is required to be strictly between 0.01 and 0.99")

  n <- nrow(data)
  numNAy <- sum(is.na(data$censor_y))
  if (numNAy>0){
    yNA.idx <- which(is.na(data$censor_y))
    data <- data[!is.na(data$censor_y),]
    message(paste("(", numNAy,
                  "observations are removed since outcome is missing)"))
  }
  regimeClass.stg1 <- as.formula(regimeClass.stg1)
  regimeClass.stg2 <- as.formula(regimeClass.stg2)
  # extract the names of the covariates in the decision rule
  p.data1 <- model.matrix(regimeClass.stg1, data)
  p.data2 <- model.matrix(regimeClass.stg2, data) # only obs.s with observed second stage treatment are contained in p.data2
  txname.stg1 <- as.character(regimeClass.stg1[[2]])
  txname.stg2 <- as.character(regimeClass.stg2[[2]])
  nvars.stg1 <- ncol(p.data1)
  nvars.stg2 <- ncol(p.data2)

  txVec1 <- try(data[, txname.stg1], silent = TRUE)
  if (is(txVec1, "try-error")) {
    stop("Variable '", paste0(txname.stg1, "' not found in 'data'."))
  }
  if (!all(unique(txVec1) %in% c(0, 1)))
    stop("The levels of treatment in the first stage should be coded as 0 or 1.")

  txVec2 <- try(data[, txname.stg2], silent = TRUE) # contain NAs
  if (is(txVec2, "try-error")) {
    stop("Variable '", paste0(txname.stg2, "' not found in 'data'."))
  }
  if (!all(unique(txVec2[which(!is.na(txVec2))]) %in% c(0, 1)))
    stop("The levels of treatment in the second stage  should be coded as 0 or 1.")

  if (moPropen1 =="BinaryRandom"){
    ph.stg1 <- rep(mean(txVec1), n)
  } else {
    moPropen1 <- as.formula(moPropen1)
    logistic.model.tx.stg1 <- glm(moPropen1, data = data, family=binomial)
    ph.stg1 <- (logistic.model.tx.stg1$fit)
  }
  if (moPropen2 =="BinaryRandom"){
    txVec2_obs <- txVec2[which(!is.na(txVec2))]
    ph.stg2 <- rep(mean(txVec2_obs), n)
    ph.stg2[which(is.na(txVec2))] <- NA
  } else { # not implemented for now.
  }

  # calculate the probability of selection
  pi.stg1 <- ph.stg1 * txVec1 + (1 - ph.stg1) * (1 - txVec1)
  pi.stg2 <- ph.stg2 * txVec2 + (1 - ph.stg2) * (1 - txVec2)

  # KM estimates of distribution of the censoring variable.
  # Assume the independent censoring assumption holds.
  # Notice that the censoring variable is missing if delta ==1.
  data$deltaC <- 1 - data$delta
  survfit_all <- survfit(Surv(censor_y, event = deltaC)~1, data=data)
  survest <- stepfun(survfit_all$time, c(1, survfit_all$surv))
  data$ghat <- survest(data$censor_y)
  g_s_Diff_Time <- survest(s_Diff_Time)
  # plot(data$censor_y, data$ghat)


  ELG <- (data$censor_y > s_Diff_Time) # the indicator of min(T(d1,null),C) <= s.
  data$w_di_vec <- rep(-999, n)
  for(i in 1:n){
    if (!ELG[i]) {
      data$w_di_vec[i] <- pi.stg1[i] * data$ghat[i]
    } else {
      data$w_di_vec[i] <- pi.stg1[i] * data$ghat[i] *  pi.stg2[i]
    }
  }

  fit <- Gene_Quantile_CenIPWE_DTR(
    data = data,
    max = max,
    tau = tau,
    s_Diff_Time = s_Diff_Time,
    regimeClass.stg1 = regimeClass.stg1,
    regimeClass.stg2 = regimeClass.stg2,
    txVec1 = txVec1,
    txVec2 = txVec2,
    nvars.stg1 = nvars.stg1,
    nvars.stg2 = nvars.stg2,
    p.data1 = p.data1,
    p.data2 = p.data2,
    Domains1 = Domains1,
    Domains2 = Domains2,
    cluster = cluster,
    pop.size = pop.size, p_level = p_level,
    s.tol = s.tol, it.num = it.num,
    sign_beta1.stg1 = sign_beta1.stg1,
    sign_beta1.stg2 = sign_beta1.stg2,
    Penalty.level = Penalty.level
  )

  fit$moPropen1 <- moPropen1
  fit$moPropen2 <- moPropen2
  fit$regimeClass.stg1 <- regimeClass.stg1
  fit$regimeClass.stg2 <- regimeClass.stg2
  fit$tau <- tau
  fit$data <- data
  fit$sign_beta1.stg1 <- sign_beta1.stg1
  fit$sign_beta1.stg2 <- sign_beta1.stg2
  fit$Penalty.level <- Penalty.level

  return(fit)
}

