#' @title Estimate the marginal quantile response of a specific dynamic TR
#'
#' @description Assume we have binary treatment options for two sequential stages with a fixed time duration between them.
#'  This means for each subject in the target population if the censored survival time or the time-to-event is beyond
#'  the timepoint of the second treatment.
#'  
#' This function evaluates a given dynamic treatment regime and returns the estimated 
#' marginal quantile response. 
#' 
#' We assume the space of two-stage treatment regimes is a cartesian product of two
#' single-stage linear treatment regime space.
#' 
#' The user facing function that applies this function is \code{\link{IPWE_Qopt_DTR_IndCen}}.
#' 
#' @inheritParams est_quant_ipwe
#' @inheritParams Gene_Quantile_CenIPWE_DTR
#' 
#' @param n the sample size
#' 
#' @param beta the vector of coefficients indexing a two-stage treatment regime
#' 
#' @param txVec2_na_omit the vector of second stage treatment for patients who indeed 
#' second stage treatment 
#' 
#' 
#' @param ELG the boolean vector of whether patients get the second stage treatment
#' 
#' @param w_di_vec the inverse probability weight for two stage experiments
#' 
#' 
#' @importFrom quantreg rq
#' @export
#' @examples
#' ##########################################################################
#' # Note: the preprocessing steps prior to calling est_quant_TwoStg_ipwe() #
#' # are wrapped up in IPWE_Qopt_DTR_IndCen().                              #
#' # w_di_vec is the inverse probability weight for two stage experiments   #
#' # We recommend users to use function IPWE_Qopt_DTR_IndCen() directly.    #
#' # Below is a simple customized calculation of the weight that only works #
#' # for this example                                                       #
#' ##########################################################################
#' 
#' library(survival)
#' # Simulate data
#' n=200
#' s_Diff_Time = 1
#' D <- simJLSDdata(n, case="a")
#' 
#' # give regime classes
#' regimeClass.stg1 <- as.formula(a0~x0)
#' regimeClass.stg2 <- as.formula(a1~x1)
#' 
#' # extract columns that matches each stage's treatment regime formula
#' p.data1 <- model.matrix(regimeClass.stg1, D)
#' 
#' # p.data2 would only contain observations with non-null value.
#' p.data2 <- model.matrix(regimeClass.stg2, D)
#' 
#' txVec1 <- D[, "a0"]
#' # get none-na second stage treatment levels in data
#' txVec2 <- D[, "a1"]
#' txVec2_na_omit <- txVec2[which(!is.na(txVec2))]
#' 
#' # Eligibility flag
#' ELG <- (D$censor_y  >  s_Diff_Time)
#' 
#' # Build weights
#' D$deltaC <- 1 - D$delta
#' survfit_all <- survfit(Surv(censor_y, event = deltaC)~1, data=D)
#' survest <- stepfun(survfit_all$time, c(1, survfit_all$surv))
#' D$ghat <- survest(D$censor_y)
#' g_s_Diff_Time <- survest(s_Diff_Time)
#' D$w_di_vec <- rep(-999, n)
#' for(i in 1:n){
#'   if (!ELG[i]) {
#'       D$w_di_vec[i] <- 0.5 * D$ghat[i]} else {
#'          D$w_di_vec[i] <- 0.5* D$ghat[i] * 0.5
#'  }
#' }
#' 
#' 
#' qhat <- est_quant_TwoStg_ipwe(n=n, beta=c(2.5,2.8), 
#'              sign_beta1.stg1 = FALSE, sign_beta1.stg2=FALSE,
#'              txVec1=txVec1, txVec2_na_omit=txVec2_na_omit, s_Diff_Time=1, 
#'              nvars.stg1=2, nvars.stg2=2, 
#'              p.data1=p.data1, 
#'              p.data2=p.data2, 
#'              censor_y=D$censor_y, 
#'              delta=D$delta, 
#'              ELG=ELG, w_di_vec=D$w_di_vec, 
#'              tau=0.3)
est_quant_TwoStg_ipwe <- function(n,
                                  beta,
                                  sign_beta1.stg1,
                                  sign_beta1.stg2,
                                  txVec1,
                                  txVec2_na_omit,
                                  s_Diff_Time,
                                  nvars.stg1,
                                  nvars.stg2,
                                  p.data1,
                                  p.data2,
                                  censor_y,
                                  delta,
                                  ELG,
                                  w_di_vec,
                                  tau,
                                  check_complete = TRUE,
                                  Penalty.level=0) {
  
  beta_stg1_full <- if (sign_beta1.stg1){ append(beta[seq(nvars.stg1 - 1)], 1, after = 1)
    } else {append(beta[seq(nvars.stg1 - 1)], -1, after = 1)}
  beta_stg2_full <- if (sign_beta1.stg2){ append(beta[- seq(nvars.stg1 - 1)], 1, after = 1)
    } else {append(beta[- seq(nvars.stg1 - 1)], -1, after = 1)}

  g1 <- (p.data1 %*% beta_stg1_full > 0)
  # p.data2 is observations with two treatments
  g2 <- (p.data2 %*% beta_stg2_full > 0)

  pre_R <- rep(TRUE, n)
  pre_R[ELG] <- (txVec2_na_omit == g2)

  R <- delta * (txVec1 == g1) * pre_R
  w_di_vec[w_di_vec < 0.01] <- 0.01
  wts <- R * (1/w_di_vec)
  fit <- rq(censor_y ~ 1, tau, weights = wts)

  if(Penalty.level==0){
    return(as.numeric(fit$coef[1]))
  } else if(Penalty.level==1) {
    return(as.numeric(fit$coef[1]) -  mean(g1))
  } else if(Penalty.level==2) {
    return(as.numeric(fit$coef[1]) -  mean(g2))
  }
}
