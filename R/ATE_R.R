#' Transporting ATE from multi-source population to an external source-specific population
#'
#' @description
#' Doubly-robust and efficient estimator for the average treatment effect of an external target population using \eqn{m} multi-source data.
#'
#' @param X The covariate matrix/data frame with \eqn{n=n_1+...+n_m} rows and q coloums.
#' @param X0 The covariate matrix/data frame with \eqn{n_0} rows and q coloums.
#' @param Y The (binary/categorical/continuous) outcome, which is a length \eqn{n} vector.
#' @param S The (numeric) source which is a length \eqn{n} vector.
#' @param A The (binary) treatment, which is a length \eqn{n} vector.
#' @param source_model The multi-nomial model for estimating \eqn{P(S=s|X)}. It has two options: \code{SL.glmnet.multinom} and \code{SL.nnet.multinom}. The default is \code{SL.glmnet.multinom}.
#' @param source_model_args The arguments (in \pkg{SuperLearner}) for the source model.
#' @param treatment_model_type The options for how the treatment_model \eqn{P(A=1|X, S=s)} is estimated. It includes \code{separate} and \code{joint}, with the default being \code{separate}. When \code{separate} is selected,
#' \eqn{P(A=1|X, S=s)} is estimated by fitting the model (regressing \eqn{A} on \eqn{X}) within each specific internal source population (S=s). When \code{joint} is selected, \eqn{P(A=1|X, S=s)}
#' is estimated by fitting the model (regressing \eqn{A} on \eqn{X} and \eqn{S}) using the multi-source population and then estimating the probability by fitting the model while suppressing the \eqn{S=s}.
#' In both cases, the propensity score is calculated as $P(A=1|X)=\sum_{s=1}^{m}P(A=1|X, S=s)P(S=s|X)$.
#' @param treatment_model The treatment model \eqn{P(A=1|X, S=s)} is estimated using \pkg{SuperLearner}. If, for example, the preference is to use only logistic regression for the probability estimation,
#' please ensure that only \code{glm} is included in the \pkg{SuperLearner} library within the \code{treatment_model_args}.
#' @param treatment_model_args The arguments (in \pkg{SuperLearner}) for the treatment model.
#' @param R_model = The R model \eqn{P(R=1|W)} is estimated using \pkg{SuperLearner}. R is a binary variable indicating the multi-source data, i.e., R is 1 if the subject belongs to the multi-source data and 0 if the subject belongs to the external data.
#' W is combination of X and X0, i.e., W=rbind(X, X0)
#' @param R_model_args = list(),
#' @param outcome_model The same as \code{treatment_model}.
#' @param outcome_model_args The arguments (in \pkg{SuperLearner}) for the outcome model.
#'
#' @details
#' Data structure: multi-source data contain outcome Y, source S, treatment A, and covariates X ($n \times p$); external data contain only covariate X0 ($n0 \times p$).
#' Once X and X0 are defined, The indicator of multi-source data, R, can be defined, i.e., R is 1 if the subject belongs to the multi-source data and 0 if the subject belongs to the external data.
#' The estimator is doubly robust and non-parametrically efficient. Three nuisance parameters are estimated,
#' the R model $q(X)=P(R=1|X)$, the propensity score model $\eta_a(X)=P(A=a|X)$, and the outcome model $\mu_a(X)=E(Y|X, A=a)$. The nuisance parameters are allowed to be estimated by \pkg{SuperLearner}. The estimator is
#' $$
#'  \dfrac{\widehat \kappa}{N}\sum\limits_{i=1}^{N} \Bigg[ I(R_i = 0) \widehat \mu_a(X_i)
#'  +I(A_i = a, R_i=1) \dfrac{1-\widehat q(X_i)}{\widehat \eta_a(X_i)\widehat q(X_i)}  \Big\{ Y_i - \widehat \mu_a(X_i) \Big\} \Bigg],
#' $$
#' where $N=n+n_0$, and $\widehat \kappa=\{N^{-1} \sum_{i=1}^N I(R_i=0)\}^{-1}$.
#' To achieve the non-parametrical efficiency and asymptotic normality, it requires that $||\widehat \mu_a(X) -\mu_a(X)||\big\{||\widehat \eta_a(X) -\eta_a(X)||+||\widehat q(X) -q(X)||\big\}=o_p(n^{-1/2})$.
#' In addition, to avoid the Donsker class assumption, the estimation is done by sample splitting and cross-fitting.
#' When one source of data is a randomized trial, it is still recommended to estimate the propensity score for optimal efficiency.
#' Since the non-parametric influence function is the same as the efficient semi-parametric efficient influence function when the propensity score is known and incorporating the assumption $Y\prep S|(X, A=a)$, the inference stays the same.
#'
#' @return An object of class "ATE_R". This object is a list with the following elements:
#'   \item{Estimate}{The point estimate of the ATE for the external data.}
#'   \item{Variance}{The asymptotic variance of the point estimate, which is calculated based on the (efficient) influence function.}
#'   \item{CI_LB}{The lower bound of the 95% confidence interval.}
#'   \item{CI_UB}{The upper bound of the 95% confidence interval.}
#'   \item{fit_outcome}{Fitted outcome model.}
#'   \item{fit_source}{Fitted source model.}
#'   \item{fit_treatment}{Fitted treatment model(s).}
#'   \item{fit_R}{Fitted R model.}
#'
#' @references Dahabreh, I.J., Robertson, S.E., Petito, L.C., Hernán, M.A. and Steingrimsson, J.A.. (2019) \emph{Efficient and robust methods for causally
#' interpretable meta‐analysis: Transporting inferences from multiple randomized trials to a target population}, Biometrics.
#'
#' @examples
#'
#' @import metafor
#' @import SuperLearner
#' @importFrom stats model.matrix predict qnorm quantile rnorm
#'
#' @export

CMetafoR.ATE.R <- function(
    X,
    X0,
    Y,
    S, # integer sequence starting from 1
    A,
    source_model = "SL.glmnet.multinom",
    source_model_args = list(),
    treatment_model_type = "separate",
    treatment_model = "SuperLearner",
    treatment_model_args = list(),
    R_model = "SuperLearner",
    R_model_args = list(),
    outcome_model = "SuperLearner",
    outcome_model_args = list()
) {
  # Total sample size
  n1 <- nrow(X)
  n0 <- nrow(X0)
  n <- n0 + n1

  # Number of sources - with format check
  unique_S <- sort(unique(S))
  no_S <- length(unique_S)
  if (! identical(seq(no_S), unique_S)) {
    stop(paste("Source", setdiff(seq(max(S)), unique_S), "is missing. "))
  }

  if (source_model %in% c("SL.glmnet.multinom", "SL.nnet.multinom")) {
    source_model_args$Y <- S
    source_model_args$X <- X
    fit_source <- do.call(what = source_model, args = source_model_args)
    PrS_X <- fit_source$pred
  } else {
    stop("Currently only support `SL.glmnet.multinom` and `SL.nnet.multinom`.")
  }

  PrA_XS <- matrix(nrow = n1, ncol = no_S)
  if (treatment_model_type == "separate") {
    fit_treatment <- vector(mode = 'list', length = no_S)
    for (s in 1:no_S) {
      id_s <- which(S == s)
      treatment_model_args$Y <- A[id_s]
      treatment_model_args$X <- X[id_s, ]
      fit_treatment_s <- do.call(what = treatment_model,
                                 args = treatment_model_args)
      PrA_XS[, s] <- predict.SuperLearner(fit_treatment_s, newdata = X)$pred
      fit_treatment[[s]] <- fit_treatment_s
    }
  } else if (treatment_model_type == "joint") {
    treatment_model_args$Y <- A
    treatment_model_args$X <- cbind(X, S)
    fit_treatment <- do.call(what = treatment_model,
                             args = treatment_model_args)
    for (s in 1:no_S) {
      PrA_XS[, s] <- predict.SuperLearner(fit_treatment,
                                          newdata = cbind(X, s))$pred
    }
  } else {
    stop("Type has to be either 'separate' or 'joint'.")
  }

  R_model_args$Y <- c(rep(1, n1), rep(0, n0))
  R_model_args$X <- rbind(X, X0)
  fit_R <- do.call(what = R_model, args = R_model_args)
  PrR_X <- predict.SuperLearner(fit_R, newdata = X)$pred

  outcome_model_args$Y <- Y
  outcome_model_args$X <- data.frame(A, X)
  fit_outcome <- do.call(what = outcome_model, args = outcome_model_args)
  pred_Y1 <- predict.SuperLearner(fit_outcome,
                                  newdata = data.frame(A = 1, X))$pred
  pred_Y0 <- predict.SuperLearner(fit_outcome,
                                  newdata = data.frame(A = 0, X))$pred
  predY_AX <- cbind(pred_Y1, pred_Y0)



  # estimators

  eta1 <- PrA_XS * PrS_X
  eta0 <- (1 - PrA_XS) * PrS_X

  pred_Y1_X0 <- predict.SuperLearner(fit_outcome,
                                     newdata = data.frame(A = 1, X0))$pred
  pred_Y0_X0 <- predict.SuperLearner(fit_outcome,
                                     newdata = data.frame(A = 0, X0))$pred

  pred_Y1_X1 <- predict.SuperLearner(fit_outcome,
                                     newdata = data.frame(A = 1, X))$pred
  pred_Y0_X1 <- predict.SuperLearner(fit_outcome,
                                     newdata = data.frame(A = 0, X))$pred

  # I_xr <- which(X0[, 1] == x_tilde)

  gamma <- n/n0 # length(I_xr)

  tmp1 <- matrix(0, nrow = n0, ncol = 2)
  tmp1[, 1] <- pred_Y1_X0   #[I_xr, ]
  tmp1[, 2] <- pred_Y0_X0   #[I_xr, ]

  tmp2 <- matrix(0, nrow = n1, ncol = 2)
  I_xa1 <- which(A == 1)
  I_xa0 <- which(A == 0)

  tmp2[I_xa1, 1] <- (1 - PrR_X[I_xa1])/PrR_X[I_xa1]/eta1[I_xa1] * (Y[I_xa1] - pred_Y1_X1[I_xa1])
  tmp2[I_xa0, 2] <- (1 - PrR_X[I_xa0])/PrR_X[I_xa0]/eta0[I_xa0] * (Y[I_xa0] - pred_Y0_X1[I_xa0])

  tmp <- colSums(rbind(tmp1, tmp2))
  phi <- gamma/n * tmp

  tmp1 <- tmp1 - rep(phi, each = n0)
  phi_var <- gamma/n^2 * colSums(rbind(tmp1, tmp2)^2)

  names(phi) <- paste0("A=", c(1, 0))

  names(phi_var) <- paste0("A=", c(1, 0))

  lb <- phi - qnorm(p = 0.975) * sqrt(phi_var)
  ub <- phi + qnorm(p = 0.975) * sqrt(phi_var)

  # tmax <- apply(abs(matrix(rnorm(length(unique(X[,1])) * 1e6),
  #                          nrow = length(unique(X[,1])), ncol = 1e6)), 2, max)
  # qtmax <- quantile(tmax, 0.95)
  #
  # lb_scb <- phi - qtmax * sqrt(phi_var)
  # ub_scb <- phi + qtmax * sqrt(phi_var)

  output <- list(Estimate = phi,
                 Variance = phi_var,
                 CI_LB = lb,
                 CI_UB = ub,
                 fit_outcome = fit_outcome,
                 fit_source = fit_source,
                 fit_treatment = fit_treatment,
                 fit_R = fit_R)
  class(output) <- 'ATE_R'

  return(output)
}


