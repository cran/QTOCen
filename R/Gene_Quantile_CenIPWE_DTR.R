#' @title A low-level function for the generic optimization step in 
#' estimating dynamic Quanilte-optimal treatment regime for censored data
#' 
#' @description This function supports wrapper functions for two stage Quantile-optimal
#' treatment regime estimation, namely
#'  \code{IPWE_Qopt_DTR_IndCen}.
#' @param data raw data.frame 
#' @param tau a quantile level of interest
#' @param max Maximization (TRUE) or Minimizing (FALSE). Determines if genoud minimizes or maximizes the objective function.
#' @param s_Diff_Time the length of time between the first stage treatment and the
#' second stage treatment
#' @param txVec1 the vector of treatment received at the first stage
#' @param txVec2 the vector of treatment received at the second stage, it expects entries
#' to be \code{NA} for patients who did not receive the second treatment
#' @param nvars.stg1 number of coeffients for the decision rule of the first stage
#' @param nvars.stg2 number of coeffients for the decision rule of the second stage
#' @param regimeClass.stg1 the class of treatment regimes for stage one
#' @param regimeClass.stg2 the class of treatment regimes for stage two
#' @param sign_beta1.stg1 Is sign of the coefficient for the first non-intercept 
#' variable for the first stage known? Default is NULL, meaning user does not have contraint on
#' the sign;
#' FALSE if the coefficient for the first continuous variable 
#' is fixed to be \code{-1}; TRUE if \code{1}. We can make the search space discrete because we employ
#' \eqn{|\beta_1| = 1} scale normalizaion.
#' 
#' @param sign_beta1.stg2 Default is NULL. Similar to \code{sign_beta1.stg1}.
#' @param p.data1 the design matrix to be used for decision in stage one
#' @param p.data2 the design matrix to be used for decision in stage two
#' @inheritParams Gene_Quantile_CenIPWE
#' @inheritParams IPWE_Qopt_DTR_IndCen
#' 
#' @import survival
#' @export
#' 
#' @examples 
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
#' txVec2 <- D[, "a1"]
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
#' \dontshow{
#' fit0  <- Gene_Quantile_CenIPWE_DTR(data=D, max=TRUE,
#'   tau=0.3,
#'   regimeClass.stg1 = regimeClass.stg1,
#'   regimeClass.stg2 = regimeClass.stg2,
#'   s_Diff_Time = s_Diff_Time,
#'   txVec1 = txVec1,
#'   txVec2 = txVec2,
#'   nvars.stg1=2,
#'   nvars.stg2=2,
#'   p.data1=p.data1,
#'   p.data2=p.data2,
#'   sign_beta1.stg1=FALSE,
#'   sign_beta1.stg2=FALSE,
#'   p_level=1,
#'   cluster=FALSE,
#'   s.tol=0.5,
#'   it.num=1,
#'   pop.size=500,
#'   Domains1 = NULL,
#'   Domains2 = NULL,
#'   Penalty.level = 0
#'   )
#'   }
#'  \donttest{
#' fit1  <- Gene_Quantile_CenIPWE_DTR(data=D, max=TRUE,
#'   tau=0.3,
#'   regimeClass.stg1 = regimeClass.stg1,
#'   regimeClass.stg2 = regimeClass.stg2,
#'   s_Diff_Time = s_Diff_Time,
#'   txVec1 = txVec1,
#'   txVec2 = txVec2,
#'   nvars.stg1=2,
#'   nvars.stg2=2,
#'   p.data1=p.data1,
#'   p.data2=p.data2,
#'   sign_beta1.stg1=FALSE,
#'   sign_beta1.stg2=NULL,
#'   p_level=1,
#'   cluster=FALSE,
#'   s.tol=1e-6,
#'   it.num=5,
#'   pop.size=6000,
#'   Domains1 = NULL,
#'   Domains2 = NULL,
#'   Penalty.level = 0
#'   )
#'  }
#' 

Gene_Quantile_CenIPWE_DTR <- function(data,
                                      max,
                                      tau,
                                      regimeClass.stg1,
                                      regimeClass.stg2,
                                      s_Diff_Time,
                                      txVec1,
                                      txVec2,
                                      nvars.stg1,
                                      nvars.stg2,
                                      p.data1,
                                      p.data2,
                                      sign_beta1.stg1,
                                      sign_beta1.stg2,
                                      p_level,
                                      cluster,
                                      s.tol,
                                      it.num,
                                      pop.size,
                                      Domains1 = NULL,
                                      Domains2 = NULL,
                                      Penalty.level = 0
)
{
  if (is.null(Domains1)) {Domains1 <- cbind(rep(-1000, nvars.stg1 - 1), rep(1000, nvars.stg1 -1))}
  if (is.null(Domains2)) {Domains2 <- cbind(rep(-1000, nvars.stg2 - 1), rep(1000, nvars.stg2 -1))}

  ELG <- (data$censor_y  >  s_Diff_Time)
  nvars <- nvars.stg1 + nvars.stg2 - 2
  est.wrapper <- function(sign_beta1.stg1, sign_beta1.stg2) {
    est <- genoud(
      fn = est_quant_TwoStg_ipwe,
      n = length(ELG),
      sign_beta1.stg1 = sign_beta1.stg1,
      sign_beta1.stg2 = sign_beta1.stg2,
      txVec1 = txVec1,
      txVec2_na_omit = txVec2[!is.na(txVec2)],
      s_Diff_Time = s_Diff_Time,
      nvars = nvars,
      nvars.stg1 = nvars.stg1,
      nvars.stg2 = nvars.stg2,
      p.data1 = p.data1,
      p.data2 = p.data2,
      Domains = rbind(Domains1, Domains2),
      censor_y = data$censor_y,
      delta = data$delta,
      ELG = ELG,
      w_di_vec = data$w_di_vec,
      Penalty.level = Penalty.level,
      tau = tau,
      print.level = p_level,
      max = max,
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

  if (!is.null(sign_beta1.stg1)) {
    if (!is.null(sign_beta1.stg2)) {
      est_c1_c2 <- est.wrapper(sign_beta1.stg1 = sign_beta1.stg1,
                         sign_beta1.stg2 = sign_beta1.stg2)
    } else {
      est_c1_n <- est.wrapper(sign_beta1.stg1 = sign_beta1.stg1,
                            sign_beta1.stg2 = FALSE)
      est_c1_p <- est.wrapper(sign_beta1.stg1 = sign_beta1.stg1,
                            sign_beta1.stg2 = TRUE)
    }
  } else{
    if (!is.null(sign_beta1.stg2)) {
      est_n_c2 <- est.wrapper(sign_beta1.stg1 = FALSE,
                            sign_beta1.stg2 = sign_beta1.stg2)
      est_p_c2 <- est.wrapper(sign_beta1.stg1 = TRUE,
                            sign_beta1.stg2 = sign_beta1.stg2)
    } else {
      est_n_n <- est.wrapper(sign_beta1.stg1 = FALSE,
                            sign_beta1.stg2 = FALSE)

      est_n_p <- est.wrapper(sign_beta1.stg1 = FALSE,
                            sign_beta1.stg2 = TRUE)

      est_p_n <- est.wrapper(sign_beta1.stg1 = TRUE,
                            sign_beta1.stg2 = FALSE)

      est_p_p <- est.wrapper(sign_beta1.stg1 = TRUE,
                            sign_beta1.stg2 = TRUE)
    }
  }

  #  set "n"= 0 to ensure it has the same mapping as that of FALSE.
  Sign.code.map <- c("n"= 0, "p" = 1, "c2" = sign_beta1.stg2, "c1" = sign_beta1.stg1)
  List_est <- grep(ls(), pattern = '^est_', value = T)
  mtx_Qhat <- NULL
  for(i in seq(length(List_est))){
    tmp_Sign_s <- strsplit(List_est[i],split = '_')[[1]][2:3]
    Sign_s <- Sign.code.map[tmp_Sign_s]
    mtx_Qhat <- rbind(mtx_Qhat, c(get(List_est[i])$value, Sign_s))
  }
  colnames(mtx_Qhat) <- c("Qhat", "sign_stg1", "sign_stg2")

  if(max){
    est <- get(List_est[which.max(mtx_Qhat[,"Qhat"])])
    sign_est_stg1 <- mtx_Qhat[which.max(mtx_Qhat[,"Qhat"]),"sign_stg1"]
    sign_est_stg2 <- mtx_Qhat[which.max(mtx_Qhat[,"Qhat"]),"sign_stg2"]
  } else {
    est <- get(List_est[which.min(mtx_Qhat[,"Qhat"])])
  }
  est$coef1 <- append(est$par[seq(nvars.stg1-1)],  (sign_est_stg1-0.5)*2, after = 1)
  est$coef2 <- append(est$par[-seq(nvars.stg1-1)], (sign_est_stg2-0.5)*2, after = 1)
  
  ######### Output the estimation ####################
  if ( max(abs(est$coef1)) == 1000)
    message("Warning: the estimated coef. for stage 1 is on the boundary of the default domain")
  if ( max(abs(est$coef2)) == 1000)
    message("Warning: the estimated coef. for stage 2 is on the boundary of the default domain")

  est$hatQ <- est$value
  est$Match_prop1 <- mean(as.numeric(p.data1 %*% est$coef1 > 0) == txVec1)
  est$Treated_prop1 <- mean(as.numeric(p.data1 %*% est$coef1 > 0) )

  return(est)
}

