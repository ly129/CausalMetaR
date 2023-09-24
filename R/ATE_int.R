#' Transporting ATE from multi-source population to an internal source-specific population
#'
#' @description
#' Doubly-robust and efficient estimator for the average treatment effects of each internal source-specific target population using \eqn{m} multi-source data.
#'
#' @param X The covariate data frame with \eqn{n=n_1+...+n_m} rows and q coloums.
#' @param Y The (binary/continuous) outcome, which is a length \eqn{n} vector.
#' @param S The (numeric) source which is a length \eqn{n} vector.
#' @param A The (binary) treatment, which is a length \eqn{n} vector.
#' @param source_model The multi-nomial model for estimating \eqn{P(S=s|X)}. It has two options: \code{glmnet.multinom} and \code{nnet.multinom}. The default is \code{glmnet.multinom}.
#' @param source_model_args The arguments (in \pkg{SuperLearner}) for the source model.
#' @param treatment_model_type The options for how the treatment_model \eqn{P(A=1|X, S=s)} is estimated. It includes \code{separate} and \code{joint}, with the default being \code{separate}. When \code{separate} is selected,
#' \eqn{P(A=1|X, S=s)} is estimated by fitting the model (regressing \eqn{A} on \eqn{X}) within each specific internal source population (S=s). When \code{joint} is selected, \eqn{P(A=1|X, S=s)}
#' is estimated by fitting the model (regressing \eqn{A} on \eqn{X} and \eqn{S}) using the multi-source population and then estimating the probability by fitting the model while suppressing the \eqn{S=s}.
#' In both cases, the propensity score is calculated as \eqn{P(A=1|X)=\sum_{s=1}^{m}P(A=1|X, S=s)P(S=s|X)}.
#' @param treatment_model The treatment model \eqn{P(A=1|X, S=s)} is estimated using \pkg{SuperLearner}. If, for example, the preference is to use only logistic regression for the probability estimation,
#' please ensure that only \code{glm} is included in the \pkg{SuperLearner} library within the \code{treatment_model_args}.
#' @param treatment_model_args The arguments (in \pkg{SuperLearner}) for the treatment model.
#' @param outcome_model The same as \code{treatment_model}.
#' @param outcome_model_args The arguments (in \pkg{SuperLearner}) for the outcome model.
#'
#' @details
#' Data structure: multi-source data contain outcome Y, source S, treatment A, and covariates X (\eqn{n \times p}).
#' The estimator is doubly robust and non-parametrically efficient. Three nuisance parameters are estimated,
#' the source model \eqn{q_s(X)=P(S=s|X)}, the propensity score model \eqn{\eta_a(X)=P(A=a|X)}, and the outcome model \eqn{\mu_a(X)=E(Y|X, A=a)}. The nuisance parameters are allowed to be estimated by \pkg{SuperLearner}. The estimator is
#' \deqn{
#'  \dfrac{\widehat \kappa}{n}\sum\limits_{i=1}^{n} \Bigg[ I(S_i = s) \widehat \mu_a(X_i)
#'  +I(A_i = a) \dfrac{\widehat q_{s}(X_i)}{\widehat \eta_a(X_i)}  \Big\{ Y_i - \widehat \mu_a(X_i) \Big\} \Bigg],
#' }
#' where \eqn{\widehat \kappa=\{n^{-1} \sum_{i=1}^n I(S_i=s)\}^{-1}}.
#' To achieve the non-parametrical efficiency and asymptotic normality, it requires that \eqn{||\widehat \mu_a(X) -\mu_a(X)||\big\{||\widehat \eta_a(X) -\eta_a(X)||+||\widehat q_s(X) -q_s(X)||\big\}=o_p(n^{-1/2})}.
#' In addition, to avoid the Donsker class assumption, the estimation is done by sample splitting and cross-fitting.
#' When one source of data is a randomized trial, it is still recommended to estimate the propensity score for optimal efficiency.
#' Since the non-parametric influence function is the same as the efficient semi-parametric efficient influence function when the propensity score is known and incorporating the assumption \eqn{Y\perp S|(X, A=a)}, the inference stays the same.
#'
#' @return An object of class "ATE_int". This object is a list with the following elements:
#'   \item{df_dif}{A data frame containing the treatment effect (mean difference) estimates for the internal populations.}
#'   \item{df_A0}{A data frame containing the potential outcome mean estimates under A = 0 for the internal populations.}
#'   \item{df_A1}{A data frame containing the potential outcome mean estimates under A = 1 for the internal populations.}
#'   \item{fit_outcome}{Fitted outcome model.}
#'   \item{fit_source}{Fitted source model.}
#'   \item{fit_treatment}{Fitted treatment model(s).}
#'
#' @references Dahabreh, I.J., Robertson, S.E., Petito, L.C., Hernán, M.A. and Steingrimsson, J.A.. (2019) \emph{Efficient and robust methods for causally
#' interpretable meta‐analysis: Transporting inferences from multiple randomized trials to a target population}, Biometrics.
#'
#' @examples
#'
#' @export

ATE_int <- function(
    X,
    Y,
    S, # integer sequence starting from 1
    A,
    source_model = "glmnet.multinom",
    source_model_args = list(),
    treatment_model_type = "separate",
    treatment_model = "SuperLearner",
    treatment_model_args = list(),
    outcome_model = "SuperLearner",
    outcome_model_args = list()
) {
  # Total sample size
  n <- nrow(X)

  # Number of sources - with format check
  unique_S <- sort(unique(S))
  no_S <- length(unique_S)

  # Error checking
  error_check(X = X, X_external = NULL, Y = Y, S = S, A = A,
              external = FALSE, ATE = TRUE)

  # Converting factor variables into dummy variables
  X <- data.frame(model.matrix(~ ., data = X)[, -1])

  if (source_model %in% c("glmnet.multinom", "nnet.multinom")) {
    source_model_args$Y <- S
    source_model_args$X <- X
    fit_source <- do.call(what = source_model, args = source_model_args)
    PrS_X <- fit_source$pred
  } else {
    stop("Currently only support `glmnet.multinom` and `nnet.multinom`.")
  }

  PrA_XS <- matrix(nrow = n, ncol = no_S)
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
    S_factor <- as.factor(S)
    S_factor <- model.matrix(~., data.frame(S_factor))[, -1]
    S_factor_names <- colnames(S_factor)

    treatment_model_args$Y <- A
    treatment_model_args$X <- data.frame(X, S_factor)
    fit_treatment <- do.call(what = treatment_model,
                             args = treatment_model_args)
    for (s in 1:no_S) {
      S_mat <- matrix(0, nrow = n, ncol = no_S - 1)
      colnames(S_mat) <- S_factor_names
      if (s > 1){
        S_mat[, s - 1] <- 1
      }
      PrA_XS[, s] <- predict.SuperLearner(fit_treatment,
                                          newdata = data.frame(X, S_mat))$pred
    }
  } else {
    stop("Type has to be either 'separate' or 'joint'.")
  }

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

  psi <- psi_var <- matrix(nrow = no_S, ncol = 2)

  for (s in unique_S) {
    tmp1 <- tmp2 <- matrix(0, nrow = n, ncol = 2)
    I_xs <- which(S == s)
    kappa <- 1/(length(I_xs)/n)
    tmp1[I_xs, 1] <- pred_Y1[I_xs]
    tmp1[I_xs, 2] <- pred_Y0[I_xs]

    I_xa1 <- which(A == 1)
    I_xa0 <- which(A == 0)

    qs <- PrS_X[, s]

    tmp2[I_xa1, 1] <- qs[I_xa1]/eta1[I_xa1] * (Y[I_xa1] - pred_Y1[I_xa1])
    tmp2[I_xa0, 2] <- qs[I_xa0]/eta0[I_xa0] * (Y[I_xa0] - pred_Y0[I_xa0])

    tmp <- tmp1 + tmp2

    psi[s, ] <- kappa * colMeans(tmp)

    tmp1[I_xs, 1] <- pred_Y1[I_xs] - psi[s, 1]
    tmp1[I_xs, 2] <- pred_Y0[I_xs] - psi[s, 2]

    psi_var[s, ] <- kappa/n^2 * colSums((tmp1 + tmp2)^2)
  }

  psi <- cbind(psi, unname(psi[, 1] - psi[, 2]))
  psi_var <- cbind(psi_var, unname(psi_var[, 1] + psi_var[, 2]))

  rownames(psi) <- rownames(psi_var) <- paste0("S = ", unique_S)
  colnames(psi) <- colnames(psi_var) <- c("A = 1", "A = 0", "Difference")

  lb <- psi - qnorm(p = 0.975) * sqrt(psi_var)
  ub <- psi + qnorm(p = 0.975) * sqrt(psi_var)

  df_dif <-
    data.frame(Source = 1:no_S,
               Estimate = psi[, 3],
               SE = psi_var[, 3],
               ci.lb = lb[, 3],
               ci.ub = ub[, 3])
  df_A0 <-
    data.frame(Source = 1:no_S,
               Estimate = psi[, 2],
               SE = psi_var[, 2],
               ci.lb = lb[, 2],
               ci.ub = ub[, 2])
  df_A1 <-
    data.frame(Source = 1:no_S,
               Estimate = psi[, 1],
               SE = psi_var[, 1],
               ci.lb = lb[, 1],
               ci.ub = ub[, 1])

  res <- list(df_dif = df_dif,
              df_A0 = df_A0,
              df_A1 = df_A1,
              fit_outcome = fit_outcome,
              fit_source = fit_source,
              fit_treatment = fit_treatment,
              outcome_model_args = outcome_model_args,
              source_model = source_model,
              treatment_model_args = treatment_model_args)
  class(res) <- 'ATE_int'

  return(res)
}


