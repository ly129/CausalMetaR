#' Transporting subgroup treatment effects (STE) from multi-source population to an internal source-specific population
#'
#' @description
#' Doubly-robust and efficient estimator for the subgroup treatments effect (STE) of each internal source-specific target population using \eqn{m} multi-source data.
#'
#' @param X The covariate data frame with \eqn{n=n_1+...+n_m} rows and q coloums. The first column of X is the categorical effect modifier (\eqn{\widetilde X}).
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
#' Data structure: multi-source data contain outcome Y, source S, treatment A, and covariates X (\eqn{n \times p}) with the first column being the categorical effect modifier (\eqn{\widetilde X}).
#' The estimator is doubly robust and non-parametrically efficient. Three nuisance parameters are estimated,
#' the source model \eqn{q_s(X)=P(S=s|X)}, the propensity score model \eqn{\eta_a(X)=P(A=a|X)}, and the outcome model \eqn{\mu_a(X)=E(Y|X, A=a)}. The nuisance parameters are allowed to be estimated by \pkg{SuperLearner}. The estimator is
#' \deqn{
#'  \dfrac{\widehat \kappa}{n}\sum\limits_{i=1}^{n} \Bigg[ I(S_i = s, \widetilde X_i=\widetilde x) \widehat \mu_a(X_i)
#'  +I(A_i = a, \widetilde X_i=\widetilde x) \dfrac{\widehat q_{s}(X_i)}{\widehat \eta_a(X_i)}  \Big\{ Y_i - \widehat \mu_a(X_i) \Big\} \Bigg],
#' }
#' where \eqn{\widehat \kappa=\{n^{-1} \sum_{i=1}^n I(S_i=s, \widetilde X_i=\widetilde x)\}^{-1}}.
#' To achieve the non-parametrical efficiency and asymptotic normality, it requires that \eqn{||\widehat \mu_a(X) -\mu_a(X)||\big\{||\widehat \eta_a(X) -\eta_a(X)||+||\widehat q_s(X) -q_s(X)||\big\}=o_p(n^{-1/2})}.
#' In addition, to avoid the Donsker class assumption, the estimation is done by sample splitting and cross-fitting.
#' When one source of data is a randomized trial, it is still recommended to estimate the propensity score for optimal efficiency.
#' Since the non-parametric influence function is the same as the efficient semi-parametric efficient influence function when the propensity score is known and incorporating the assumption \eqn{Y\perp S|(X, A=a)}, the inference stays the same.
#'
#' @return An object of class "STE_int". This object is a list with the following elements:
#'   \item{df_dif}{A data frame containing the subgroup treatment effect (mean difference) estimates for the internal populations.}
#'   \item{df_A0}{A data frame containing the subgroup potential outcome mean estimates under A = 0 for the internal populations.}
#'   \item{df_A1}{A data frame containing the subgroup potential outcome mean estimates under A = 1 for the internal populations.}
#'   \item{fit_outcome}{Fitted outcome model.}
#'   \item{fit_source}{Fitted source model.}
#'   \item{fit_treatment}{Fitted treatment model(s).}
#'   \item{...}{Some additional elements.}
#'
#' @examples
#'
#' @export
STE_int_cf <- function(
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
    outcome_model_args = list(),
    replications = 10L
) {
  # Total sample size
  n <- nrow(X)

  # Number of sources
  unique_S <- sort(unique(S))
  no_S <- length(unique_S)

  # Error checking
  error_check(X = X, X_external = NULL, Y = Y, S = S, A = A,
              external = FALSE, ATE = FALSE)

  # Converting factor variables into dummy variables
  X1_name <- names(X)[1]
  X1 <- X[, 1]
  X <- data.frame(model.matrix(~ ., data = X)[, -1])
  unique_X <- sort(unique(X1))
  no_x_tilde <- length(unique_X)

  ## sample splitting and cross fitting loop
  K <- 4L
  phi_array <- phi_var_array <- array(dim = c(no_S, 3, no_x_tilde, K, replications))
  for (r in 1:replications) {
    ### assign k in 0, 1, 2, 3 to each individual
    id_by_SX <- partition <- vector(mode = "list", length = no_S)
    for (s in unique_S) {
      for (i in unique_X) {
        idsx <- which(X1 == i & S == s)
        nsx <- length(idsx)
        partition[[s]][[i]] <- sample(rep(seq(K) - 1, length.out = nsx))
        id_by_SX[[s]][[i]] <- idsx
      }
    }
    for (k in 1:K) {
      test.id <- sm.id <- tm.id <- om.id <- integer()
      for (s in unique_S) {
        id_S <- id_by_SX[[s]]
        part_S <- partition[[s]]
        for (i in unique_X) {
          id_SX <- id_S[[i]]
          part_SX <- part_S[[i]]
          test.id <- c(test.id, id_SX[which(part_SX == k %% K)])
          sm.id <- c(sm.id, id_SX[which(part_SX == (k + 1) %% K)])
          tm.id <- c(tm.id, id_SX[which(part_SX == (k + 2) %% K)])
          om.id <- c(om.id, id_SX[which(part_SX == (k + 3) %% K)])
        }
      }
      X_test <- X[test.id, ]
      Y_test <- Y[test.id]
      A_test <- A[test.id]
      S_test <- S[test.id]
      X1_test <- X1[test.id]

      X_sm <- X[sm.id, ]
      Y_sm <- Y[sm.id]
      A_sm <- A[sm.id]
      S_sm <- S[sm.id]

      X_tm <- X[tm.id, ]
      Y_tm <- Y[tm.id]
      A_tm <- A[tm.id]
      S_tm <- S[tm.id]

      X_om <- X[om.id, ]
      Y_om <- Y[om.id]
      A_om <- A[om.id]
      S_om <- S[om.id]

      ## source model
      if (source_model %in% c("glmnet.multinom", "nnet.multinom")) {
        source_model_args$Y <- S_sm
        source_model_args$X <- X_sm
        source_model_args$newX <- X_test
        fit_source <- do.call(what = source_model, args = source_model_args)
        PrS_X <- fit_source$pred
      } else {
        stop("Currently only support `glmnet.multinom` and `nnet.multinom`.")
      }

      ## treatment model
      PrA_XS <- matrix(nrow = nrow(X_test), ncol = no_S)
      if (treatment_model_type == "separate") {
        fit_treatment <- vector(mode = 'list', length = no_S)
        for (s in 1:no_S) {
          id_s <- which(S_tm == s)
          treatment_model_args$Y <- A_tm[id_s]
          treatment_model_args$X <- X_tm[id_s, ]
          fit_treatment_s <- do.call(what = treatment_model,
                                     args = treatment_model_args)
          PrA_XS[, s] <- predict.SuperLearner(fit_treatment_s, newdata = X_test)$pred
          fit_treatment[[s]] <- fit_treatment_s
        }
      } else if (treatment_model_type == "joint") {
        S_factor <- as.factor(S_tm)
        S_factor <- model.matrix(~., data.frame(S_factor))[, -1]
        S_factor_names <- colnames(S_factor)

        treatment_model_args$Y <- A_tm
        treatment_model_args$X <- data.frame(X_tm, S_factor)
        fit_treatment <- do.call(what = treatment_model,
                                 args = treatment_model_args)
        for (s in 1:no_S) {
          S_mat <- matrix(0, nrow = nrow(X_test), ncol = no_S - 1)
          colnames(S_mat) <- S_factor_names
          if (s > 1){
            S_mat[, s - 1] <- 1
          }
          PrA_XS[, s] <- predict.SuperLearner(fit_treatment,
                                              newdata = data.frame(X_test, S_mat))$pred
        }
      } else {
        stop("Type has to be either 'separate' or 'joint'.")
      }

      ## outcome model
      outcome_model_args$Y <- Y_om
      outcome_model_args$X <- data.frame(A_om, X_om)
      fit_outcome <- do.call(what = outcome_model, args = outcome_model_args)
      pred_Y1 <- predict.SuperLearner(fit_outcome,
                                      newdata = data.frame(A_om = 1, X_test))$pred
      pred_Y0 <- predict.SuperLearner(fit_outcome,
                                      newdata = data.frame(A_om = 0, X_test))$pred
      predY_AX <- cbind(pred_Y1, pred_Y0)

      # estimators
      eta1 <- rowSums(PrA_XS * PrS_X)
      eta0 <- rowSums((1 - PrA_XS) * PrS_X)


      output <- vector(mode = "list", length = no_x_tilde)
      names(output) <- paste(X1_name, "=", unique_X)

      for (i in 1:no_x_tilde) {
        phi <- phi_var <- matrix(nrow = no_S, ncol = 2)
        x_tilde <- unique_X[i]
        for (s in unique_S) {
          tmp1 <- tmp2 <- matrix(0, nrow = nrow(X_test), ncol = 2)
          I_xs <- which((X1_test == x_tilde) & (S_test == s))
          kappa <- 1/(length(I_xs)/nrow(X_test))
          tmp1[I_xs, 1] <- pred_Y1[I_xs]
          tmp1[I_xs, 2] <- pred_Y0[I_xs]

          I_xa1 <- which((X1_test == x_tilde) & (A_test == 1))
          I_xa0 <- which((X1_test == x_tilde) & (A_test == 0))

          qs <- PrS_X[, s]

          tmp2[I_xa1, 1] <- qs[I_xa1]/eta1[I_xa1] * (Y_test[I_xa1] - pred_Y1[I_xa1])
          tmp2[I_xa0, 2] <- qs[I_xa0]/eta0[I_xa0] * (Y_test[I_xa0] - pred_Y0[I_xa0])

          tmp <- tmp1 + tmp2

          phi[s, ] <- kappa * colMeans(tmp)

          if (is.infinite(phi[s, 1])) stop("r", r, ", k", k, ", xtilde", i)

          tmp1[I_xs, 1] <- pred_Y1[I_xs] - phi[s, 1]
          tmp1[I_xs, 2] <- pred_Y0[I_xs] - phi[s, 2]

          phi_var[s, ] <- kappa/nrow(X_test)^2 * colSums((tmp1 + tmp2)^2)
        }
        phi <- cbind(phi, unname(phi[, 1] - phi[, 2]))
        phi_var <- cbind(phi_var, unname(phi_var[, 1] + phi_var[, 2]))

        phi_array[, , i, k, r] <- phi
        phi_var_array[, , i, k, r] <- phi_var
      }
    }
  }

  phi_cf <- apply(apply(phi_array, MARGIN = c(1, 2, 3, 5), FUN = mean), MARGIN = 1:3, FUN = median)
  phi_var_cf <- apply(apply(phi_var_array, MARGIN = c(1, 2, 3, 5), FUN = mean), MARGIN = 1:3, FUN = median)

  qt <- qnorm(p = 0.975)
  qtmax <- quantile(apply(abs(matrix(rnorm(no_x_tilde * 1e6),
                                     nrow = no_x_tilde, ncol = 1e6)), 2, max), 0.95)

  df_A0 <- df_A1 <- df_dif <- matrix(nrow = no_S * no_x_tilde, ncol = 8L)
  df_A0[, 1] <- df_A1[, 1] <- df_dif[, 1] <- rep(unique_S, each = no_x_tilde)
  df_A0[, 2] <- df_A1[, 2] <- df_dif[, 2] <- rep(unique_X, no_S)

  for (i in 1:(no_S * no_x_tilde)) {
    S_id <- df_A0[i, 1]
    X_id <- df_A0[i, 2]
    # Estimates
    df_A0[i, 3] <- phi_cf[S_id, 2, X_id]
    df_A1[i, 3] <- phi_cf[S_id, 1, X_id]
    df_dif[i, 3] <- phi_cf[S_id, 3, X_id]
    # SE
    df_A0[i, 4] <- sqrt(phi_var_cf[S_id, 2, X_id])
    df_A1[i, 4] <- sqrt(phi_var_cf[S_id, 1, X_id])
    df_dif[i, 4] <- sqrt(phi_var_cf[S_id, 3, X_id])
  }
  # lb
  df_A0[, 5] <- df_A0[, 3] - qt * df_A0[, 4]
  df_A1[, 5] <- df_A1[, 3] - qt * df_A1[, 4]
  df_dif[, 5] <- df_dif[, 3] - qt * df_dif[, 4]
  # ub
  df_A0[, 6] <- df_A0[, 3] + qt * df_A0[, 4]
  df_A1[, 6] <- df_A1[, 3] + qt * df_A1[, 4]
  df_dif[, 6] <- df_dif[, 3] + qt * df_dif[, 4]
  # lb_scb
  df_A0[, 7] <- df_A0[, 3] - qtmax * df_A0[, 4]
  df_A1[, 7] <- df_A1[, 3] - qtmax * df_A1[, 4]
  df_dif[, 7] <- df_dif[, 3] - qtmax * df_dif[, 4]
  # ub_scb
  df_A0[, 8] <- df_A0[, 3] + qtmax * df_A0[, 4]
  df_A1[, 8] <- df_A1[, 3] + qtmax * df_A1[, 4]
  df_dif[, 8] <- df_dif[, 3] + qtmax * df_dif[, 4]

  df_A0 <- as.data.frame(df_A0)
  df_A1 <- as.data.frame(df_A1)
  df_dif <- as.data.frame(df_dif)

  names(df_A0) <- names(df_A1) <- names(df_dif) <- c("Source",
                                                     "Subgroup",
                                                     "Estimate",
                                                     "SE",
                                                     "ci.lb",
                                                     "ci.ub",
                                                     "scb.lb",
                                                     "scb.ub")

  snames <- character(length = no_x_tilde * no_S)
  snames[1:(no_x_tilde * no_S) %% no_x_tilde == 1] <- c(paste("Source", unique_S))
  xtildenames <- rep(paste(X1_name, "=", unique_X), no_S)

  res <- list(
    df_dif = df_dif,
    df_A0 = df_A0,
    df_A1 = df_A1,
    outcome_model_args = outcome_model_args,
    source_model = source_model,
    treatment_model_args = treatment_model_args
  )

  res$snames <- snames
  res$xtildenames <- xtildenames

  class(res) <- 'STE_int'

  return(res)
}


