#' Transporting ATE from multi-source population to an internal source-specific population
#'
#' @description
#' Doubly-robust and efficient estimator for the average treatment effects of each internal source-specific target population using \eqn{m} multi-source data.
#'
#' @param X The covariate data frame with \eqn{n=n_1+...+n_{|S|}} rows and \eqn{p} columns. Character variables will be converted to factors. Numeric variables will be used as is.
#' @param Y The length \eqn{n} outcome vector.
#' @param S The source indicator which is a length \eqn{n} vector or factor. If \code{S} is a factor, it will maintain its level order, otherwise it will be converted to a factor with default level order. The order will be carried over to the outputs and plots.
#' @param A The binary treatment (1 for treated and 0 for untreated), which is a length \eqn{n} vector.
#' @param cross_fitting Logical, indicating whether sample splitting and cross fitting procedure should be used.
#' @param replications Integer, the number of sample splitting and cross fitting replications to performe, if \code{cross_fitting = TRUE}. Default is \code{10L}.
#' @param source_model The (penalized) multinomial logistic regression for estimating \eqn{P(S=s|X)}. It has two options: "\code{MN.glmnet}" (default) and "\code{MN.nnet}", which use \pkg{glmnet} and \pkg{nnet} respectively.
#' @param source_model_args The arguments (in \pkg{SuperLearner}) for the source model.
#' @param treatment_model_type How the propensity score \eqn{P(A=1|X)=\sum_{s \in S} P(A=1|X, S=s)P(S=s|X)} is estimated. Options include "\code{separate}" (default) and "\code{joint}". If "\code{separate}", \eqn{P(A=1|X, S=s)} is estimated by regressing \eqn{A} on \eqn{X} within each specific internal source population \eqn{S=s}. If "\code{joint}", \eqn{P(A=1|X, S=s)} is estimated by regressing \eqn{A} on \eqn{X} and \eqn{S} using the multi-source population.
#' @param treatment_model_args The arguments (in \pkg{SuperLearner}) for the treatment model.
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
#' @return An object of class "ATE_nested". This object is a list with the following elements:
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

ATE_nested <- function(
    X, # predictor matrix
    Y, # outcome
    S, # source
    A, # treatment
    cross_fitting = FALSE,
    replications = 10L,
    source_model = "MN.glmnet",
    source_model_args = list(),
    treatment_model_type = "separate",
    treatment_model_args = list(),
    outcome_model_args = list()
) {
  # For future possibilities
  treatment_model <- outcome_model <- "SuperLearner"

  # # Error checking
  # error_check(X = X, X_external = NULL, Y = Y, S = S, A = A,
  #             external = FALSE, ATE = FALSE)

  # Total sample size
  n <- nrow(X)

  # Number of sources
  if (!is.factor(S)) S <- factor(S)
  unique_S <- levels(S)
  no_S <- length(unique_S)
  # S with one-hot encoding to prevent potential problems in SuperLearner
  S_dummy <- data.frame(model.matrix(~ 0 + S)[, -1])
  colnames(S_dummy) <- unique_S[-1]

  # Converting character variables into factor variables
  X <- sapply(X, FUN = function(xx) {
    if (is.character(xx)) xx <- factor(xx)
  })
  # Converting factor variables into dummy variables
  X <- data.frame(model.matrix(~ ., data = X)[, -1])

  if (cross_fitting) {
    ## sample splitting and cross fitting loop
    K <- 4L
    phi_array_cf <- phi_se_array_cf <- array(dim = c(no_S, 3, K, replications),
                                             dimnames = list(unique_S,
                                                             c("A = 1", "A = 0", "Difference"),
                                                             NULL,
                                                             NULL))
    for (r in 1:replications) {
      ### assign k in 0, 1, 2, 3 to each individual
      id_by_S <- partition <- vector(mode = "list", length = no_S)
      for (s in unique_S) {
        ids <- which(S == s)
        ns <- length(ids)
        partition[[s]] <- sample(rep(seq(K) - 1, length.out = ns))
        id_by_S[[s]] <- ids
      }
      for (k in 1:K) {
        test_id <- sm_id <- tm_id <- om_id <- integer()
        for (s in unique_S) {
          id_S <- id_by_S[[s]]
          part_S <- partition[[s]]
          test_id <- c(test_id, id_S[part_S == k %% K])
          sm_id <- c(sm_id, id_S[part_S == (k + 1) %% K])
          tm_id <- c(tm_id, id_S[part_S == (k + 2) %% K])
          om_id <- c(om_id, id_S[part_S == (k + 3) %% K])
        }
        X_test <- X[test_id, ]
        Y_test <- Y[test_id]
        A_test <- A[test_id]
        S_test <- S[test_id]

        X_sm <- X[sm_id, ]
        S_sm <- S[sm_id]

        X_tm <- X[tm_id, ]
        A_tm <- A[tm_id]
        S_tm <- S[tm_id]
        if (treatment_model_type == "joint") S_dummy_tm <- S_dummy[tm_id, ]

        X_om <- X[om_id, ]
        Y_om <- Y[om_id]
        A_om <- A[om_id]

        ## source model
        if (source_model %in% c("MN.glmnet", "MN.nnet")) {
          source_model_args$Y <- S_sm
          source_model_args$X <- X_sm
          source_model_args$newX <- X_test
          fit_source <- do.call(what = source_model, args = source_model_args)
          PrS_X <- fit_source$pred
        } else {
          stop("Currently only support `MN.glmnet` and `MN.nnet`.")
        }

        ## treatment model
        PrA_XS <- matrix(nrow = length(test_id), ncol = no_S)
        colnames(PrA_XS) <- unique_S
        if (treatment_model_type == "separate") {
          fit_treatment <- vector(mode = 'list', length = no_S)
          for (s in unique_S) {
            id_s <- which(S_tm == s)
            treatment_model_args$Y <- A_tm[id_s]
            treatment_model_args$X <- X_tm[id_s, ]
            fit_treatment_s <- do.call(what = treatment_model,
                                       args = treatment_model_args)
            PrA_XS[, s] <- predict.SuperLearner(fit_treatment_s, newdata = X_test)$pred
            fit_treatment[[s]] <- fit_treatment_s
          }
        } else if (treatment_model_type == "joint") {
          treatment_model_args$Y <- A_tm
          treatment_model_args$X <- data.frame(S_dummy_tm, X_tm)
          fit_treatment <- do.call(what = treatment_model,
                                   args = treatment_model_args)
          for (s in unique_S) {
            S_mat <- matrix(0, nrow = length(test_id), ncol = no_S - 1)
            colnames(S_mat) <- unique_S[-1]
            if (s %in% unique_S[-1]) S_mat[, s] <- 1
            PrA_XS[, s] <- predict.SuperLearner(fit_treatment,
                                                newdata = data.frame(S_mat, X_test))$pred
          }
        } else {
          stop("The 'treatment_model_type' has to be either 'separate' or 'joint'.")
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

        ## estimators
        eta1 <- rowSums(PrA_XS * PrS_X)
        eta0 <- rowSums((1 - PrA_XS) * PrS_X)

        phi <- phi_var <- matrix(nrow = no_S, ncol = 2)
        rownames(phi) <- rownames(phi_var) <- unique_S
        for (s in unique_S) {
          tmp1 <- tmp2 <- matrix(0, nrow = nrow(X_test), ncol = 2)
          I_xs <- which(S_test == s)
          kappa <- 1/(length(I_xs)/nrow(X_test))
          tmp1[I_xs, 1] <- pred_Y1[I_xs]
          tmp1[I_xs, 2] <- pred_Y0[I_xs]

          I_xa1 <- which(A_test == 1)
          I_xa0 <- which(A_test == 0)

          qs <- PrS_X[, s]

          tmp2[I_xa1, 1] <- qs[I_xa1]/eta1[I_xa1] * (Y_test[I_xa1] - pred_Y1[I_xa1])
          tmp2[I_xa0, 2] <- qs[I_xa0]/eta0[I_xa0] * (Y_test[I_xa0] - pred_Y0[I_xa0])

          tmp <- tmp1 + tmp2

          phi[s, ] <- kappa * colMeans(tmp)

          tmp1[I_xs, 1] <- pred_Y1[I_xs] - phi[s, 1]
          tmp1[I_xs, 2] <- pred_Y0[I_xs] - phi[s, 2]

          phi_var[s, ] <- kappa/nrow(X_test)^2 * colSums((tmp1 + tmp2)^2)
        }

        phi_array_cf[, , k, r] <- cbind(phi, phi[, 1] - phi[, 2])
        phi_se_array_cf[, , k, r] <- sqrt(cbind(phi_var, phi_var[, 1] + phi_var[, 2]))
      }
    }

    phi_array <- apply(apply(phi_array_cf,
                             MARGIN = c(1, 2, 4),
                             FUN = mean),
                       MARGIN = 1:2,
                       FUN = median)
    phi_se_array <- apply(apply(phi_se_array_cf,
                                MARGIN = c(1, 2, 4),
                                FUN = mean),
                          MARGIN = 1:2,
                          FUN = median)
    # end of ATE_nested with cross-fitting
  } else {
    # start of regular ATE_nested
    if (source_model %in% c("MN.glmnet", "MN.nnet")) {
      source_model_args$Y <- S
      source_model_args$X <- X
      fit_source <- do.call(what = source_model, args = source_model_args)
      PrS_X <- fit_source$pred
    } else {
      stop("Currently only support `MN.glmnet` and `MN.nnet`.")
    }

    PrA_XS <- matrix(nrow = n, ncol = no_S)
    colnames(PrA_XS) <- unique_S
    if (treatment_model_type == "separate") {
      fit_treatment <- vector(mode = 'list', length = no_S)
      for (s in unique_S) {
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
      treatment_model_args$X <- data.frame(S_dummy, X)
      fit_treatment <- do.call(what = treatment_model,
                               args = treatment_model_args)
      for (s in unique_S) {
        S_mat <- matrix(0, nrow = n, ncol = no_S - 1)
        colnames(S_mat) <- unique_S[-1]
        if (s %in% unique_S[-1]) S_mat[, s] <- 1
        PrA_XS[, s] <- predict.SuperLearner(fit_treatment,
                                            newdata = data.frame(S_mat, X))$pred
      }
    } else {
      stop("The 'treatment_model_type' has to be either 'separate' or 'joint'.")
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
    eta1 <- rowSums(PrA_XS * PrS_X)
    eta0 <- rowSums((1 - PrA_XS) * PrS_X)

    phi_array <- phi_se_array <- array(dim = c(no_S, 3),
                                       dimnames = list(unique_S,
                                                       c("A = 1", "A = 0", "Difference")))
    phi <- phi_var <- matrix(nrow = no_S, ncol = 2)
    rownames(phi) <- rownames(phi_var) <- unique_S
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

      phi[s, ] <- kappa * colMeans(tmp)

      tmp1[I_xs, 1] <- pred_Y1[I_xs] - phi[s, 1]
      tmp1[I_xs, 2] <- pred_Y0[I_xs] - phi[s, 2]

      phi_var[s, ] <- kappa/n^2 * colSums((tmp1 + tmp2)^2)
    }

    phi_array[, ] <- cbind(phi, phi[, 1] - phi[, 2])
    phi_se_array[, ] <- sqrt(cbind(phi_var, phi_var[, 1] + phi_var[, 2]))
  }

  lb <- phi_array - qnorm(p = 0.975) * phi_se_array
  ub <- phi_array + qnorm(p = 0.975) * phi_se_array

  df_dif <- df_A0 <- df_A1 <-
    data.frame(Source = unique_S,
               Estimate = NA,
               SE = NA,
               ci.lb = NA,
               ci.ub = NA)

  for (s in unique_S) {
    rowid <- which(df_dif$Source == s)
    df_dif$Estimate[rowid] <- phi_array[s, "Difference"]
    df_dif$SE[rowid] <- phi_se_array[s, "Difference"]
    df_dif$ci.lb[rowid] <- lb[s, "Difference"]
    df_dif$ci.ub[rowid] <- ub[s, "Difference"]

    df_A1$Estimate[rowid] <- phi_array[s, "A = 1"]
    df_A1$SE[rowid] <- phi_se_array[s, "A = 1"]
    df_A1$ci.lb[rowid] <- lb[s, "A = 1"]
    df_A1$ci.ub[rowid] <- ub[s, "A = 1"]

    df_A0$Estimate[rowid] <- phi_array[s, "A = 0"]
    df_A0$SE[rowid] <- phi_se_array[s, "A = 0"]
    df_A0$ci.lb[rowid] <- lb[s, "A = 0"]
    df_A0$ci.ub[rowid] <- ub[s, "A = 0"]
  }

  if (cross_fitting) {
    fit_outcome <- fit_source <- fit_treatment <- NULL
  }

  res <- list(df_dif = df_dif,
              df_A0 = df_A0,
              df_A1 = df_A1,
              fit_outcome = fit_outcome,
              fit_source = fit_source,
              fit_treatment = fit_treatment,
              outcome_model_args = outcome_model_args,
              source_model = source_model,
              treatment_model_args = treatment_model_args)

  res$source_names <- unique_S

  class(res) <- 'ATE_nested'

  return(res)
}


