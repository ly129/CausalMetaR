#' Transporting subgroup treatment effects (STE) from multi-source population to an external source-specific population
#'
#' @description
#' Doubly-robust and efficient estimator for the subgroup treatment effects (STE) of an external target population using \eqn{m} multi-source data.
#'
#' @param X The covariate matrix/data frame with \eqn{n=n_1+...+n_m} rows and q coloums. The first column of X is the categorical effect modifier (\eqn{\widetilde X}).
#' @param X0 The covariate matrix/data frame with \eqn{n_0} rows and q coloums. The first column of X is the categorical effect modifier (\eqn{\widetilde X}).
#' @param Y The (binary/categorical/continuous) outcome, which is a length \eqn{n} vector.
#' @param S The (numeric) source which is a length \eqn{n} vector.
#' @param A The (binary) treatment, which is a length \eqn{n} vector.
#' @param EM The effect modifier (FILL IN DETAILS)
#' @param source_model The multi-nomial model for estimating \eqn{P(S=s|X)}. It has two options: \code{SL.glmnet.multinom} and \code{SL.nnet.multinom}. The default is \code{SL.glmnet.multinom}.
#' @param source_model_args The arguments (in \pkg{SuperLearner}) for the source model.
#' @param treatment_model_type The options for how the treatment_model \eqn{P(A=1|X, S=s)} is estimated. It includes \code{separate} and \code{joint}, with the default being \code{separate}. When \code{separate} is selected,
#' \eqn{P(A=1|X, S=s)} is estimated by fitting the model (regressing \eqn{A} on \eqn{X}) within each specific internal source population (S=s). When \code{joint} is selected, \eqn{P(A=1|X, S=s)}
#' is estimated by fitting the model (regressing \eqn{A} on \eqn{X} and \eqn{S}) using the multi-source population and then estimating the probability by fitting the model while suppressing the \eqn{S=s}.
#' In both cases, the propensity score is calculated as \eqn{P(A=1|X)=\sum_{s=1}^{m}P(A=1|X, S=s)P(S=s|X)}.
#' @param treatment_model The treatment model \eqn{P(A=1|X, S=s)} is estimated using \pkg{SuperLearner}. If, for example, the preference is to use only logistic regression for the probability estimation,
#' please ensure that only \code{glm} is included in the \pkg{SuperLearner} library within the \code{treatment_model_args}.
#' @param treatment_model_args The arguments (in \pkg{SuperLearner}) for the treatment model.
#' @param external_model The R model \eqn{P(R=1|W)} is estimated using \pkg{SuperLearner}. R is a binary variable indicating the multi-source data, i.e., R is 1 if the subject belongs to the multi-source data and 0 if the subject belongs to the external data.
#' W is combination of X and X0, i.e., W=rbind(X, X0)
#' @param external_model_args = list(),
#' @param outcome_model The same as \code{treatment_model}.
#' @param outcome_model_args The arguments (in \pkg{SuperLearner}) for the outcome model.
#'
#' @details
#' Data structure: multi-source data contain outcome Y, source S, treatment A, and covariates X (\eqn{n \times p}); external data contain only covariate X0 (\eqn{n_0 \times p}).
#' Once X and X0 are defined, The indicator of multi-source data, R, can be defined, i.e., R is 1 if the subject belongs to the multi-source data and 0 if the subject belongs to the external data.
#' The estimator is doubly robust and non-parametrically efficient. Three nuisance parameters are estimated,
#' the R model \eqn{q(X)=P(R=1|X)}, the propensity score model \eqn{\eta_a(X)=P(A=a|X)}, and the outcome model \eqn{\mu_a(X)=E(Y|X, A=a)}. The nuisance parameters are allowed to be estimated by \pkg{SuperLearner}. The estimator is
#' \deqn{
#'  \dfrac{\widehat \kappa}{N}\sum\limits_{i=1}^{N} \Bigg[ I(R_i = 0) \widehat \mu_a(X_i)
#'  +I(A_i = a, R_i=1) \dfrac{1-\widehat q(X_i)}{\widehat \eta_a(X_i)\widehat q(X_i)}  \Big\{ Y_i - \widehat \mu_a(X_i) \Big\} \Bigg],
#' }
#' where \eqn{N=n+n_0}, and \eqn{\widehat \kappa=\{N^{-1} \sum_{i=1}^N I(R_i=0)\}^{-1}}.
#' To achieve the non-parametrical efficiency and asymptotic normality, it requires that \eqn{||\widehat \mu_a(X) -\mu_a(X)||\big\{||\widehat \eta_a(X) -\eta_a(X)||+||\widehat q(X) -q(X)||\big\}=o_p(n^{-1/2})}.
#' In addition, to avoid the Donsker class assumption, the estimation is done by sample splitting and cross-fitting.
#' When one source of data is a randomized trial, it is still recommended to estimate the propensity score for optimal efficiency.
#' Since the non-parametric influence function is the same as the efficient semi-parametric efficient influence function when the propensity score is known and incorporating the assumption \eqn{Y\perp S|(X, A=a)}, the inference stays the same.
#'
#' @return An object of class "STE_ext". This object is a list with the following elements:
#'   \item{Estimates}{The point estimates of the STE for the external data.}
#'   \item{Variances}{The asymptotic variances of the point estimates, which is calculated based on the (efficient) influence function.}
#'   \item{CI_LB}{The lower bounds of the 95% confidence intervals.}
#'   \item{CI_UB}{The upper bounds of the 95% confidence intervals.}
#'   \item{SCB_LB}{The lower bounds of the 95% simultaneous confidence bands.}
#'   \item{SCB_UB}{The upper bounds of the 95% simultaneous confidence bands.}
#'   \item{fit_outcome}{Fitted outcome model.}
#'   \item{fit_source}{Fitted source model.}
#'   \item{fit_treatment}{Fitted treatment model(s).}
#'   \item{fit_external}{Fitted external model.}
#'
#' @examples
#'
#' @export

STE_ext <- function(
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
    external_model = "SuperLearner",
    external_model_args = list(),
    outcome_model = "SuperLearner",
    outcome_model_args = list(),
    EM
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
    treatment_model_args$X <- data.frame(X, S)
    fit_treatment <- do.call(what = treatment_model,
                             args = treatment_model_args)
    for (s in 1:no_S) {
      PrA_XS[, s] <- predict.SuperLearner(fit_treatment,
                                          newdata = data.frame(X, S = s))$pred
    }
  } else {
    stop("Type has to be either 'separate' or 'joint'.")
  }

  external_model_args$Y <- c(rep(1, n1), rep(0, n0))
  external_model_args$X <- rbind(X, X0)
  fit_external <- do.call(what = external_model, args = external_model_args)
  PrR_X <- predict.SuperLearner(fit_external, newdata = X)$pred

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

  I_xr <- which(X0[, 1] == EM)

  gamma <- n/length(I_xr)

  tmp1 <- matrix(0, nrow = n0, ncol = 2)
  tmp1[I_xr, 1] <- pred_Y1_X0[I_xr, ]
  tmp1[I_xr, 2] <- pred_Y0_X0[I_xr, ]

  tmp2 <- matrix(0, nrow = n1, ncol = 2)
  I_xa1 <- which((X[, 1] == EM) & (A == 1))
  I_xa0 <- which((X[, 1] == EM) & (A == 0))

  tmp2[I_xa1, 1] <- (1 - PrR_X[I_xa1])/PrR_X[I_xa1]/eta1[I_xa1] * (Y[I_xa1] - pred_Y1_X1[I_xa1])
  tmp2[I_xa0, 2] <- (1 - PrR_X[I_xa0])/PrR_X[I_xa0]/eta0[I_xa0] * (Y[I_xa0] - pred_Y0_X1[I_xa0])

  tmp <- colSums(rbind(tmp1, tmp2))
  phi <- gamma/n * tmp

  tmp1[I_xr, ] <- tmp1[I_xr, ] - rep(phi, each = length(I_xr))
  phi_var <- gamma/n^2 * colSums(rbind(tmp1, tmp2)^2)

  names(phi) <- paste0("A = ", c(1, 0))

  names(phi_var) <- paste0("A = ", c(1, 0))

  lb <- phi - qnorm(p = 0.975) * sqrt(phi_var)
  ub <- phi + qnorm(p = 0.975) * sqrt(phi_var)

  tmax <- apply(abs(matrix(rnorm(length(unique(X[,1])) * 1e6),
                           nrow = length(unique(X[,1])), ncol = 1e6)), 2, max)
  qtmax <- quantile(tmax, 0.95)

  lb_scb <- phi - qtmax * sqrt(phi_var)
  ub_scb <- phi + qtmax * sqrt(phi_var)

  # STE
  plot_phi <- unname(phi[1] - phi[2])
  plot_phi_var <- unname(phi_var[1] + phi_var[2])
  plot_phi_CI <- plot_phi + c(-1, 1) * qnorm(p = 0.975) * sqrt(plot_phi_var)
  plot_phi_SCB <- plot_phi + c(-1, 1) * qtmax * sqrt(plot_phi_var)

  output <- list(Estimates = phi,
                 Variances = phi_var,
                 CI_LB = lb,
                 CI_UB = ub,
                 SCB_LB = lb_scb,
                 SCB_UB = ub_scb,
                 fit_outcome = fit_outcome,
                 fit_source = fit_source,
                 fit_treatment = fit_treatment,
                 fit_external = fit_external,
                 plot_phi = plot_phi,
                 plot_phi_var = plot_phi_var,
                 plot_phi_CI = plot_phi_CI,
                 plot_phi_SCB = plot_phi_SCB)
  class(output) <- 'STE_ext'

  return(output)
}


