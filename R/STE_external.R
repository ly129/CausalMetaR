#' Transporting subgroup treatment effects (STE) from multi-source population to an external source-specific population
#'
#' @description
#' Doubly-robust and efficient estimator for the subgroup treatment effects (STE) of an external target population using \eqn{m} multi-source data.
#'
#' @param X The covariate data frame with \eqn{n=n_1+...+n_m} rows and q coloums. The first column of X is the categorical effect modifier (\eqn{\widetilde X}).
#' @param X_external The covariate matrix/data frame with \eqn{n_0} rows and q coloums. The first column of X is the categorical effect modifier (\eqn{\widetilde X}).
#' @param Y The (binary/continuous) outcome, which is a length \eqn{n} vector.
#' @param S The (numeric) source which is a length \eqn{n} vector.
#' @param A The (binary) treatment, which is a length \eqn{n} vector.
#' @param source_model The multi-nomial model for estimating \eqn{P(S=s|X)}. It has two options: \code{MN.glmnet} and \code{MN.nnet}. The default is \code{MN.glmnet}.
#' @param source_model_args The arguments (in \pkg{SuperLearner}) for the source model.
#' @param treatment_model_type The options for how the treatment_model \eqn{P(A=1|X, S=s)} is estimated. It includes \code{separate} and \code{joint}, with the default being \code{separate}. When \code{separate} is selected,
#' \eqn{P(A=1|X, S=s)} is estimated by fitting the model (regressing \eqn{A} on \eqn{X}) within each specific internal source population (S=s). When \code{joint} is selected, \eqn{P(A=1|X, S=s)}
#' is estimated by fitting the model (regressing \eqn{A} on \eqn{X} and \eqn{S}) using the multi-source population and then estimating the probability by fitting the model while suppressing the \eqn{S=s}.
#' In both cases, the propensity score is calculated as \eqn{P(A=1|X)=\sum_{s=1}^{m}P(A=1|X, S=s)P(S=s|X)}.
#' @param treatment_model The treatment model \eqn{P(A=1|X, S=s)} is estimated using \pkg{SuperLearner}. If, for example, the preference is to use only logistic regression for the probability estimation,
#' please ensure that only \code{glm} is included in the \pkg{SuperLearner} library within the \code{treatment_model_args}.
#' @param treatment_model_args The arguments (in \pkg{SuperLearner}) for the treatment model.
#' @param external_model The R model \eqn{P(R=1|W)} is estimated using \pkg{SuperLearner}. R is a binary variable indicating the multi-source data, i.e., R is 1 if the subject belongs to the multi-source data and 0 if the subject belongs to the external data.
#' W is combination of X and X_external, i.e., W=rbind(X, X_external)
#' @param external_model_args = list(),
#' @param outcome_model The same as \code{treatment_model}.
#' @param outcome_model_args The arguments (in \pkg{SuperLearner}) for the outcome model.
#'
#' @details
#' Data structure: multi-source data contain outcome Y, source S, treatment A, and covariates X (\eqn{n \times p}); external data contain only covariate X_external (\eqn{n_0 \times p}).
#' Once X and X_external are defined, The indicator of multi-source data, R, can be defined, i.e., R is 1 if the subject belongs to the multi-source data and 0 if the subject belongs to the external data.
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
#' @return An object of class "STE_external". This object is a list with the following elements:
#'   \item{df_dif}{A data frame containing the subgroup treatment effect (mean difference) estimates for the extenal data.}
#'   \item{df_A0}{A data frame containing the subgroup potential outcome mean estimates under A = 0 for the extenal data.}
#'   \item{df_A1}{A data frame containing the subgroup potential outcome mean estimates under A = 1 for the extenal data.}
#'   \item{fit_outcome}{Fitted outcome model.}
#'   \item{fit_source}{Fitted source model.}
#'   \item{fit_treatment}{Fitted treatment model(s).}
#'   \item{fit_external}{Fitted external model.}
#'
#' @examples
#'
#' @export

STE_external <- function(
    X,
    X_external,
    EM,
    EM_external,
    Y,
    S,
    A,
    cross_fitting = FALSE,
    replications = 10L,
    source_model = "MN.glmnet",
    source_model_args = list(),
    treatment_model_type = "separate",
    treatment_model_args = list(),
    external_model_args = list(),
    outcome_model_args = list()
) {
  # For future possibilities
  treatment_model <- external_model <- outcome_model <- "SuperLearner"

  # # Error checking
  # error_check(X = X, X_external = X_external, Y = Y, S = S, A = A,
  #             external = TRUE, ATE = FALSE)

  # Total sample size
  n1 <- nrow(X)
  n0 <- nrow(X_external)
  n <- n0 + n1

  # Number of sources
  if (!is.factor(S)) S <- factor(S)
  unique_S <- levels(S)
  no_S <- length(unique_S)
  # S with one-hot encoding to prevent potential problems in SuperLearner
  S_dummy <- data.frame(model.matrix(~ 0 + S)[, -1])
  colnames(S_dummy) <- unique_S[-1]

  # Number of EM
  if (!is.factor(EM)) EM <- factor(EM)
  unique_EM <- levels(EM)
  no_EM <- length(unique_EM)
  # EM with one-hot encoding to prevent potential problems in SuperLearner
  EM_dummy <- data.frame(model.matrix(~ 0 + EM)[, -1])
  colnames(EM_dummy) <- unique_EM[-1]

  # Number of EM_external
  if (!is.factor(EM_external)) EM_external <- factor(EM_external)
  unique_EM_ext <- levels(EM_external)
  no_EM_ext <- length(unique_EM_ext)
  # EM_external with one-hot encoding to prevent potential problems in SuperLearner
  EM_ext_dummy <- data.frame(model.matrix(~ 0 + EM_external)[, -1])
  colnames(EM_ext_dummy) <- unique_EM_ext[-1]

  if (no_EM != no_EM_ext) stop("Different number of subgroups in external data.")

  # Converting character variables into factor variables
  X <- sapply(X, FUN = function(xx) {
    if (is.character(xx)) xx <- factor(xx)
  })
  X_external <- sapply(X_external, FUN = function(xx) {
    if (is.character(xx)) xx <- factor(xx)
  })

  # Converting factor variables into dummy variables
  X <- data.frame(model.matrix(~ ., data = X)[, -1])
  X_external <- data.frame(model.matrix(~ ., data = X_external)[, -1])

  if (cross_fitting) {
    K <- 5L
    psi_array_cf <- psi_se_array_cf <- array(dim = c(no_EM, 3, K, replications),
                                             dimnames = list(unique_EM,
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
      partition_ext <- sample(rep(1:K, length.out = n0))
      for (k in 1:K) {
        test_id <- sm_id <- tm_id <- om_id <- xm_id <- integer()
        for (s in unique_S) {
          id_S <- id_by_S[[s]]
          part_S <- partition[[s]]
          test_id <- c(test_id, id_S[part_S == k %% K])
          sm_id <- c(sm_id, id_S[part_S == (k + 1) %% K])
          tm_id <- c(tm_id, id_S[part_S == (k + 2) %% K])
          om_id <- c(om_id, id_S[part_S == (k + 3) %% K])
          xm_id <- c(xm_id, id_S[part_S == (k + 4) %% K])
        }
        ext_xm_id <- which(partition_ext == k)

        X_test <- X[test_id, ]
        Y_test <- Y[test_id]
        A_test <- A[test_id]
        EM_dummy_test <- EM_dummy[test_id, ]
        EM_test <- EM[test_id]

        X_sm <- X[sm_id, ]
        S_sm <- S[sm_id]
        EM_dummy_sm <- EM_dummy[sm_id, ]

        X_tm <- X[tm_id, ]
        A_tm <- A[tm_id]
        S_tm <- S[tm_id]
        if (treatment_model_type == "joint") S_dummy_tm <- S_dummy[tm_id, ]
        EM_dummy_tm <- EM_dummy[tm_id, ]

        X_om <- X[om_id, ]
        Y_om <- Y[om_id]
        A_om <- A[om_id]
        EM_dummy_om <- EM_dummy[om_id, ]

        X_xm <- X[xm_id, ]
        EM_dummy_xm <- EM_dummy[xm_id, ]
        X_ext_xm <- X_external[ext_xm_id, ]
        EM_ext_dummy_xm <- EM_ext_dummy[ext_xm_id, ]

        # source model
        if (source_model %in% c("MN.glmnet", "MN.nnet")) {
          source_model_args$Y <- S_sm
          source_model_args$X <- data.frame(EM_dummy_sm, X_sm)
          source_model_args$newX <- data.frame(EM_dummy_test, X_test)
          fit_source <- do.call(what = source_model, args = source_model_args)
          PrS_X <- fit_source$pred
        } else {
          stop("Currently only support `MN.glmnet` and `MN.nnet`.")
        }

        # treatment model
        PrA_XS <- matrix(nrow = length(test_id), ncol = no_S)
        colnames(PrA_XS) <- unique_S
        if (treatment_model_type == "separate") {
          fit_treatment <- vector(mode = 'list', length = no_S)
          for (s in unique_S) {
            id_s <- which(S_tm == s)
            treatment_model_args$Y <- A_tm[id_s]
            treatment_model_args$X <- data.frame(EM_dummy_tm, X_tm)[id_s, ]
            fit_treatment_s <- do.call(what = treatment_model,
                                       args = treatment_model_args)
            PrA_XS[, s] <- predict.SuperLearner(fit_treatment_s,
                                                newdata = data.frame(EM_dummy_test, X_test))$pred
            fit_treatment[[s]] <- fit_treatment_s
          }
        } else if (treatment_model_type == "joint") {
          treatment_model_args$Y <- A_tm
          treatment_model_args$X <- data.frame(S_dummy_tm, EM_dummy_tm, X_tm)
          fit_treatment <- do.call(what = treatment_model,
                                   args = treatment_model_args)
          for (s in unique_S) {
            S_mat <- matrix(0, nrow = length(test_id), ncol = no_S - 1)
            colnames(S_mat) <- unique_S[-1]
            if (s %in% unique_S[-1]) S_mat[, s] <- 1
            PrA_XS[, s] <- predict.SuperLearner(fit_treatment,
                                                newdata = data.frame(S_mat, EM_dummy_test, X_test))$pred
          }
        } else {
          stop("The 'treatment_model_type' has to be either 'separate' or 'joint'.")
        }

        # external model
        external_model_args$Y <- c(rep(1, length(xm_id)), rep(0, length(ext_xm_id)))
        external_model_args$X <- rbind(data.frame(EM_dummy_xm, X_xm),
                                       data.frame(EM_ext_dummy_xm, X_ext_xm))
        fit_external <- do.call(what = external_model, args = external_model_args)
        PrR_X <- predict.SuperLearner(fit_external, newdata = data.frame(EM_dummy_test, X_test))$pred

        # outcome model
        outcome_model_args$Y <- Y_om
        outcome_model_args$X <- data.frame(A = A_om, EM_dummy_om, X_om)
        fit_outcome <- do.call(what = outcome_model, args = outcome_model_args)
        pred_Y1 <- predict.SuperLearner(fit_outcome,
                                        newdata = data.frame(A = 1, EM_dummy_test, X_test))$pred
        pred_Y0 <- predict.SuperLearner(fit_outcome,
                                        newdata = data.frame(A = 0, EM_dummy_test, X_test))$pred
        predY_AX <- cbind(pred_Y1, pred_Y0)

        # estimators
        eta1 <- rowSums(PrA_XS * PrS_X)
        eta0 <- rowSums((1 - PrA_XS) * PrS_X)

        pred_Y1_X0 <- c(predict.SuperLearner(fit_outcome,
                                             newdata = data.frame(A = 1, EM_ext_dummy_xm, X_ext_xm))$pred)
        pred_Y0_X0 <- c(predict.SuperLearner(fit_outcome,
                                             newdata = data.frame(A = 0, EM_ext_dummy_xm, X_ext_xm))$pred)
        pred_Y1_X1 <- c(predict.SuperLearner(fit_outcome,
                                             newdata = data.frame(A = 1, EM_dummy_test, X_test))$pred)
        pred_Y0_X1 <- c(predict.SuperLearner(fit_outcome,
                                             newdata = data.frame(A = 0, EM_dummy_test, X_test))$pred)

        n_kr <- length(ext_xm_id) + length(test_id)
        for (em in unique_EM) {
          I_xr <- which(EM_external[ext_xm_id] == em)

          gamma <- n_kr/length(I_xr)

          tmp1 <- matrix(0, nrow = length(ext_xm_id), ncol = 2)
          tmp1[I_xr, 1] <- pred_Y1_X0[I_xr]
          tmp1[I_xr, 2] <- pred_Y0_X0[I_xr]

          tmp2 <- matrix(0, nrow = length(test_id), ncol = 2)
          I_xa1 <- which((EM_test == em) & (A_test == 1))
          I_xa0 <- which((EM_test == em) & (A_test == 0))

          tmp2[I_xa1, 1] <- (1 - PrR_X[I_xa1])/PrR_X[I_xa1]/eta1[I_xa1] * (Y_test[I_xa1] - pred_Y1_X1[I_xa1])
          tmp2[I_xa0, 2] <- (1 - PrR_X[I_xa0])/PrR_X[I_xa0]/eta0[I_xa0] * (Y_test[I_xa0] - pred_Y0_X1[I_xa0])

          tmp <- colSums(rbind(tmp1, tmp2))
          psi <- gamma/n_kr * tmp

          tmp1[I_xr, ] <- tmp1[I_xr, ] - rep(psi, each = length(I_xr))
          psi_var <- gamma/n_kr^2 * colSums(rbind(tmp1, tmp2)^2)

          psi_array_cf[em, , k, r] <- c(psi, psi[1] - psi[2])
          psi_se_array_cf[em, , k, r] <- sqrt(c(psi_var, psi_var[2] + psi_var[2]))
        }
      }
    }

    psi_array <- apply(apply(psi_array_cf,
                             MARGIN = c(1, 2, 4),
                             FUN = mean),
                       MARGIN = 1:2,
                       FUN = median)
    psi_se_array <- apply(apply(psi_se_array_cf,
                                MARGIN = c(1, 2, 4),
                                FUN = mean),
                          MARGIN = 1:2,
                          FUN = median)
    # end of STE_external with cross-fitting
  } else {
    # start of regular STE_external
    if (source_model %in% c("MN.glmnet", "MN.nnet")) {
      source_model_args$Y <- S
      source_model_args$X <- data.frame(EM_dummy, X)
      fit_source <- do.call(what = source_model, args = source_model_args)
      PrS_X <- fit_source$pred
    } else {
      stop("Currently only support `MN.glmnet` and `MN.nnet`.")
    }

    PrA_XS <- matrix(nrow = n1, ncol = no_S)
    colnames(PrA_XS) <- unique_S
    if (treatment_model_type == "separate") {
      fit_treatment <- vector(mode = 'list', length = no_S)
      for (s in unique_S) {
        id_s <- which(S == s)
        treatment_model_args$Y <- A[id_s]
        treatment_model_args$X <- data.frame(EM_dummy, X)[id_s, ]
        fit_treatment_s <- do.call(what = treatment_model,
                                   args = treatment_model_args)
        PrA_XS[, s] <- predict.SuperLearner(fit_treatment_s,
                                            newdata = data.frame(EM_dummy, X))$pred
        fit_treatment[[s]] <- fit_treatment_s
      }
    } else if (treatment_model_type == "joint") {
      treatment_model_args$Y <- A
      treatment_model_args$X <- data.frame(S_dummy, EM_dummy, X)
      fit_treatment <- do.call(what = treatment_model,
                               args = treatment_model_args)
      for (s in unique_S) {
        S_mat <- matrix(0, nrow = n1, ncol = no_S - 1)
        colnames(S_mat) <- unique_S[-1]
        if (s %in% unique_S[-1]) S_mat[, s] <- 1
        PrA_XS[, s] <- predict.SuperLearner(fit_treatment,
                                            newdata = data.frame(S_mat, EM_dummy, X))$pred
      }
    } else {
      stop("The 'treatment_model_type' has to be either 'separate' or 'joint'.")
    }

    external_model_args$Y <- c(rep(1, n1), rep(0, n0))
    external_model_args$X <- rbind(data.frame(EM_dummy, X),
                                   data.frame(EM_ext_dummy, X_external))
    fit_external <- do.call(what = external_model, args = external_model_args)
    PrR_X <- predict.SuperLearner(fit_external, newdata = data.frame(EM_dummy, X))$pred

    outcome_model_args$Y <- Y
    outcome_model_args$X <- data.frame(A, EM_dummy, X)
    fit_outcome <- do.call(what = outcome_model, args = outcome_model_args)
    pred_Y1 <- predict.SuperLearner(fit_outcome,
                                    newdata = data.frame(A = 1, EM_dummy, X))$pred
    pred_Y0 <- predict.SuperLearner(fit_outcome,
                                    newdata = data.frame(A = 0, EM_dummy, X))$pred
    predY_AX <- cbind(pred_Y1, pred_Y0)

    # estimators
    eta1 <- rowSums(PrA_XS * PrS_X)
    eta0 <- rowSums((1 - PrA_XS) * PrS_X)

    pred_Y1_X0 <- predict.SuperLearner(fit_outcome,
                                       newdata = data.frame(A = 1, EM_ext_dummy, X_external))$pred
    pred_Y0_X0 <- predict.SuperLearner(fit_outcome,
                                       newdata = data.frame(A = 0, EM_ext_dummy, X_external))$pred

    pred_Y1_X1 <- predict.SuperLearner(fit_outcome,
                                       newdata = data.frame(A = 1, EM_dummy, X))$pred
    pred_Y0_X1 <- predict.SuperLearner(fit_outcome,
                                       newdata = data.frame(A = 0, EM_dummy, X))$pred

    psi_array <- psi_se_array <- array(dim = c(no_EM, 3),
                                       dimnames = list(unique_EM,
                                                       c("A = 1", "A = 0", "Difference")))

    psi <- psi_var <- matrix(nrow = no_EM, ncol = 2)
    rownames(psi) <- rownames(psi_var) <- unique_EM
    for (em in unique_EM) {
      I_xr <- which(EM_external == em)

      gamma <- n/length(I_xr)

      tmp1 <- matrix(0, nrow = n0, ncol = 2)
      tmp1[I_xr, 1] <- pred_Y1_X0[I_xr, ]
      tmp1[I_xr, 2] <- pred_Y0_X0[I_xr, ]

      tmp2 <- matrix(0, nrow = n1, ncol = 2)
      I_xa1 <- which((EM == em) & (A == 1))
      I_xa0 <- which((EM == em) & (A == 0))

      tmp2[I_xa1, 1] <- (1 - PrR_X[I_xa1])/PrR_X[I_xa1]/eta1[I_xa1] * (Y[I_xa1] - pred_Y1_X1[I_xa1])
      tmp2[I_xa0, 2] <- (1 - PrR_X[I_xa0])/PrR_X[I_xa0]/eta0[I_xa0] * (Y[I_xa0] - pred_Y0_X1[I_xa0])

      tmp <- colSums(rbind(tmp1, tmp2))
      psi[em, ] <- gamma/n * tmp

      tmp1[I_xr, ] <- tmp1[I_xr, ] - rep(psi[em, ], each = length(I_xr))
      psi_var[em, ] <- gamma/n^2 * colSums(rbind(tmp1, tmp2)^2)
    }

    psi_array[, ] <- cbind(psi, psi[, 1] - psi[, 2])
    psi_se_array[, ] <- sqrt(cbind(psi_var, psi_var[, 1] + psi_var[, 2]))
  }

  qt <- qnorm(p = 0.975)
  qtmax <- quantile(apply(abs(matrix(rnorm(no_EM * 1e6),
                                     nrow = no_EM, ncol = 1e6)), 2, max), 0.95)

  lb <- psi_array - qnorm(p = 0.975) * sqrt(psi_se_array)
  ub <- psi_array + qnorm(p = 0.975) * sqrt(psi_se_array)

  lb_scb <- psi_array - qtmax * sqrt(psi_se_array)
  ub_scb <- psi_array + qtmax * sqrt(psi_se_array)

  # Put results in a data frame
  df_dif <- df_A1 <- df_A0 <-
    data.frame(Subgroup = unique_EM,
               Estimate = NA,
               SE = NA,
               ci.lb = NA,
               ci.ub = NA,
               scb.lb = NA,
               scb.ub = NA)

  for (em in unique_EM) {
    rowid <- which(df_dif$Subgroup == em)
    df_dif$Estimate[rowid] <- psi_array[em, "Difference"]
    df_dif$SE[rowid] <- psi_se_array[em, "Difference"]
    df_dif$ci.lb[rowid] <- lb[em, "Difference"]
    df_dif$ci.ub[rowid] <- ub[em, "Difference"]
    df_dif$scb.lb[rowid] <- lb_scb[em, "Difference"]
    df_dif$scb.ub[rowid] <- ub_scb[em, "Difference"]

    df_A1$Estimate[rowid] <- psi_array[em, "A = 1"]
    df_A1$SE[rowid] <- psi_se_array[em, "A = 1"]
    df_A1$ci.lb[rowid] <- lb[em, "A = 1"]
    df_A1$ci.ub[rowid] <- ub[em, "A = 1"]
    df_A1$scb.lb[rowid] <- lb_scb[em, "A = 1"]
    df_A1$scb.ub[rowid] <- ub_scb[em, "A = 1"]

    df_A0$Estimate[rowid] <- psi_array[em, "A = 0"]
    df_A0$SE[rowid] <- psi_se_array[em, "A = 0"]
    df_A0$ci.lb[rowid] <- lb[em, "A = 0"]
    df_A0$ci.ub[rowid] <- ub[em, "A = 0"]
    df_A0$scb.lb[rowid] <- lb_scb[em, "A = 0"]
    df_A0$scb.ub[rowid] <- ub_scb[em, "A = 0"]
  }

  if (cross_fitting) {
    fit_outcome <- fit_source <- fit_treatment <- fit_external <- NULL
  }

  res <- list(
    df_dif = df_dif,
    df_A0 = df_A0,
    df_A1 = df_A1,
    fit_outcome = fit_outcome,
    fit_source = fit_source,
    fit_treatment = fit_treatment,
    fit_external = fit_external,
    outcome_model_args = outcome_model_args,
    source_model = source_model,
    treatment_model_args = treatment_model_args,
    external_model_args = external_model_args
  )

  class(res) <- 'STE_external'

  return(res)
}


