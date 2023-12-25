#' Transporting subgroup treatment effects (STE) from multi-source population to a nested population
#'
#' @description
#' Doubly-robust and efficient estimator for the subgroup treatments effect (STE) of each nested target population using \eqn{m} multi-source data.
#'
#' @param X The covariate data frame with \eqn{n=n_1+...+n_{|S|}} rows and \eqn{p} columns. Character variables will be converted to factors. Numeric variables will be used as is.
#' @param EM The effect modifier, which is a vector or factor of length \eqn{n}. If \code{EM} is a factor, it will maintain its subgroup level order, otherwise it will be converted to a factor with default level order.
#' @param Y The length \eqn{n} outcome vector.
#' @param S The source indicator which is a length \eqn{n} vector or factor. If \code{S} is a factor, it will maintain its level order, otherwise it will be converted to a factor with default level order. The order will be carried over to the outputs and plots.
#' @param A The binary treatment (1 for treated and 0 for untreated), which is a length \eqn{n} vector.
#' @param cross_fitting Logical, indicating whether sample splitting and cross fitting procedure should be used.
#' @param replications Integer, the number of sample splitting and cross fitting replications to perform, if \code{cross_fitting = TRUE}. Default is \code{10L}.
#' @param source_model The (penalized) multinomial logistic regression for estimating \eqn{P(S=s|X)}. It has two options: "\code{MN.glmnet}" (default) and "\code{MN.nnet}", which use \pkg{glmnet} and \pkg{nnet} respectively.
#' @param source_model_args The arguments (in \pkg{glmnet} or \pkg{nnet}) for the source model.
#' @param treatment_model_type How the propensity score \eqn{P(A=1|X)=\sum_{s \in S} P(A=1|X, S=s)P(S=s|X)} is estimated. Options include "\code{separate}" (default) and "\code{joint}". If "\code{separate}", \eqn{P(A=1|X, S=s)} is estimated by regressing \eqn{A} on \eqn{X} within each specific nested population \eqn{S=s}. If "\code{joint}", \eqn{P(A=1|X, S=s)} is estimated by regressing \eqn{A} on \eqn{X} and \eqn{S} using the multi-source population.
#' @param treatment_model_args The arguments (in \pkg{SuperLearner}) for the treatment model.
#' @param outcome_model_args The arguments (in \pkg{SuperLearner}) for the outcome model.
#' @param show_progress Logical, indicating whether to print a progress bar for the cross-fit replicates completed, if \code{cross_fitting = TRUE}.
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
#' To achieve non-parametric efficiency and asymptotic normality, it requires that \eqn{||\widehat \mu_a(X) -\mu_a(X)||\big\{||\widehat \eta_a(X) -\eta_a(X)||+||\widehat q_s(X) -q_s(X)||\big\}=o_p(n^{-1/2})}.
#' In addition, to avoid the Donsker class assumption, the estimation is done by sample splitting and cross-fitting.
#' When one source of data is a randomized trial, it is still recommended to estimate the propensity score for optimal efficiency.
#' Since the non-parametric influence function is the same as the efficient semi-parametric efficient influence function when the propensity score is known and incorporating the assumption \eqn{Y\perp S|(X, A=a)}, the inference stays the same.
#'
#' @return An object of class "STE_nested". This object is a list with the following elements:
#'   \item{df_dif}{A data frame containing the subgroup treatment effect (mean difference) estimates for the nested populations.}
#'   \item{df_A0}{A data frame containing the subgroup potential outcome mean estimates under A = 0 for the nested populations.}
#'   \item{df_A1}{A data frame containing the subgroup potential outcome mean estimates under A = 1 for the nested populations.}
#'   \item{fit_outcome}{Fitted outcome model.}
#'   \item{fit_source}{Fitted source model.}
#'   \item{fit_treatment}{Fitted treatment model(s).}
#'   \item{...}{Some additional elements.}
#'
#' @examples
#' sn <- STE_nested(
#'   X = dat_nested[, 2:10],
#'   Y = dat_nested$Y,
#'   EM = dat_nested$EM,
#'   S = dat_nested$S,
#'   A = dat_nested$A,
#'   cross_fitting = FALSE,
#'   source_model = "MN.nnet",
#'   source_model_args = list(),
#'   treatment_model_type = "separate",
#'   treatment_model_args = list(
#'     family = binomial(),
#'     SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
#'     cvControl = list(V = 5L)
#'   ),
#'   outcome_model_args = list(
#'     family = gaussian(),
#'     SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
#'     cvControl = list(V = 5L)
#'   )
#' )
#' @export
STE_nested <- function(
    X, # predictor matrix
    Y, # outcome
    EM, # effect modifier
    S, # source
    A, # treatment
    cross_fitting = FALSE,
    replications = 10L,
    source_model = "MN.glmnet",
    source_model_args = list(),
    treatment_model_type = "separate",
    treatment_model_args = list(),
    outcome_model_args = list(),
    show_progress = TRUE
) {
  # For future possibilities
  treatment_model <- outcome_model <- "SuperLearner"

  # Total sample size
  n <- nrow(X)

  # Checking lengths of all data inputs
  if (length(EM) != n |length(Y) != n | length(S) != n | length(A) != n){
    stop(paste0('All data inputs must have the same length:\n',
                'X has length ', n, '.\n',
                'EM has length ', length(EM), '.\n',
                'Y has length ', length(Y), '.\n',
                'S has length ', length(S), '.\n',
                'A has length ', length(A), '.'))
  }

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

  # Converting character variables into factor variables
  X <- as.data.frame(sapply(X, FUN = function(xx) {
    if (is.character(xx)) {
      factor(xx)
    } else {
      xx
    }
  }))
  # Converting factor variables into dummy variables
  X <- data.frame(model.matrix(~ ., data = X)[, -1])

  # Checking treatment data
  if (!is.numeric(A)){
    stop('A must be a numeric vector')
  }
  if (!all(A %in% c(0, 1))){
    stop('A must only take values 0 or 1')
  }

  if (source_model == "MN.nnet" & !('trace' %in% source_model_args)){
    source_model_args$trace <- FALSE
  }

  if (cross_fitting) {
    if (show_progress){
      pb <- progress::progress_bar$new(total = replications,
                                       clear = FALSE,
                                       format = 'Cross-fitting progress [:bar] :percent, Elapsed time :elapsed, Est. time remaining :eta')
    }
    ## sample splitting and cross fitting loop
    K <- 4L
    phi_array_cf <- phi_se_array_cf <- array(dim = c(no_S, 3, no_EM, K, replications),
                                             dimnames = list(unique_S,
                                                             c("A = 1", "A = 0", "Difference"),
                                                             unique_EM,
                                                             NULL,
                                                             NULL))
    for (r in 1:replications) {
      ### assign k in 0, 1, 2, 3 to each individual
      id_by_SX <- partition <- vector(mode = "list", length = no_S)
      for (s in unique_S) {
        for (em in unique_EM) {
          idsx <- which(EM == em & S == s)
          nsx <- length(idsx)
          partition[[s]][[em]] <- sample(rep(seq(K) - 1, length.out = nsx))
          id_by_SX[[s]][[em]] <- idsx
        }
      }
      for (k in 1:K) {
        test_id <- sm_id <- tm_id <- om_id <- integer()
        for (s in unique_S) {
          id_S <- id_by_SX[[s]]
          part_S <- partition[[s]]
          for (em in unique_EM) {
            id_SX <- id_S[[em]]
            part_SX <- part_S[[em]]
            test_id <- c(test_id, id_SX[part_SX == k %% K])
            sm_id <- c(sm_id, id_SX[part_SX == (k + 1) %% K])
            tm_id <- c(tm_id, id_SX[part_SX == (k + 2) %% K])
            om_id <- c(om_id, id_SX[part_SX == (k + 3) %% K])
          }
        }
        X_test <- X[test_id, ]
        Y_test <- Y[test_id]
        A_test <- A[test_id]
        S_test <- S[test_id]
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

        ## source model
        if (source_model %in% c("MN.glmnet", "MN.nnet")) {
          source_model_args$Y <- S_sm
          source_model_args$X <- data.frame(EM_dummy_sm, X_sm)
          source_model_args$newX <- data.frame(EM_dummy_test, X_test)
          fit_source <- do.call(what = source_model, args = source_model_args)
          PrS_X <- fit_source$pred
        } else {
          stop("Currently only support `MN.glmnet` and `MN.nnet`.")
        }

        ## treatment model
        PrA_XS <- matrix(nrow = nrow(X_test), ncol = no_S)
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
            S_mat <- matrix(0, nrow = nrow(X_test), ncol = no_S - 1)
            colnames(S_mat) <- unique_S[-1]
            if (s %in% unique_S[-1]) S_mat[, s] <- 1
            PrA_XS[, s] <- predict.SuperLearner(fit_treatment,
                                                newdata = data.frame(S_mat, EM_dummy_test, X_test))$pred
          }
        } else {
          stop("The 'treatment_model_type' has to be either 'separate' or 'joint'.")
        }

        ## outcome model
        outcome_model_args$Y <- Y_om
        outcome_model_args$X <- data.frame(A_om, EM_dummy_om, X_om)
        fit_outcome <- do.call(what = outcome_model, args = outcome_model_args)
        pred_Y1 <- predict.SuperLearner(fit_outcome,
                                        newdata = data.frame(A_om = 1, EM_dummy_test, X_test))$pred
        pred_Y0 <- predict.SuperLearner(fit_outcome,
                                        newdata = data.frame(A_om = 0, EM_dummy_test, X_test))$pred
        predY_AX <- cbind(pred_Y1, pred_Y0)

        # estimators
        eta1 <- rowSums(PrA_XS * PrS_X)
        eta0 <- rowSums((1 - PrA_XS) * PrS_X)

        for (em in unique_EM) {
          phi <- phi_var <- matrix(nrow = no_S, ncol = 2)
          rownames(phi) <- rownames(phi_var) <- unique_S
          for (s in unique_S) {
            tmp1 <- tmp2 <- matrix(0, nrow = nrow(X_test), ncol = 2)
            I_xs <- which((EM_test == em) & (S_test == s))
            kappa <- 1/(length(I_xs)/nrow(X_test))
            tmp1[I_xs, 1] <- pred_Y1[I_xs]
            tmp1[I_xs, 2] <- pred_Y0[I_xs]

            I_xa1 <- which((EM_test == em) & (A_test == 1))
            I_xa0 <- which((EM_test == em) & (A_test == 0))

            qs <- PrS_X[, s]

            tmp2[I_xa1, 1] <- qs[I_xa1]/eta1[I_xa1] * (Y_test[I_xa1] - pred_Y1[I_xa1])
            tmp2[I_xa0, 2] <- qs[I_xa0]/eta0[I_xa0] * (Y_test[I_xa0] - pred_Y0[I_xa0])

            tmp <- tmp1 + tmp2

            phi[s, ] <- kappa * colMeans(tmp)

            # if (is.infinite(phi[s, 1])) stop("r", r, ", k", k, ", xtilde", i)

            tmp1[I_xs, 1] <- pred_Y1[I_xs] - phi[s, 1]
            tmp1[I_xs, 2] <- pred_Y0[I_xs] - phi[s, 2]

            phi_var[s, ] <- kappa/nrow(X_test)^2 * colSums((tmp1 + tmp2)^2)
          }

          phi_array_cf[, , em, k, r] <- cbind(phi, phi[, 1] - phi[, 2])
          phi_se_array_cf[, , em, k, r] <- sqrt(cbind(phi_var, phi_var[, 1] + phi_var[, 2]))
        }
      }
      if (show_progress){
        pb$tick()
      }
    }

    phi_array <- apply(apply(phi_array_cf,
                             MARGIN = c(1, 2, 3, 5),
                             FUN = mean),
                       MARGIN = 1:3,
                       FUN = median)
    phi_se_array <- apply(apply(phi_se_array_cf,
                                MARGIN = c(1, 2, 3, 5),
                                FUN = mean),
                          MARGIN = 1:3,
                          FUN = median)
    # end of STE_nested with cross-fitting
  } else {
    # start of regular STE_nested
    if (source_model %in% c("MN.glmnet", "MN.nnet")) {
      source_model_args$Y <- S
      source_model_args$X <- data.frame(EM_dummy, X)
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
        treatment_model_args$X <- data.frame(EM_dummy, X)[id_s, ]
        fit_treatment_s <- do.call(what = treatment_model,
                                   args = treatment_model_args)
        PrA_XS[, s] <- predict.SuperLearner(fit_treatment_s, newdata = data.frame(EM_dummy, X))$pred
        fit_treatment[[s]] <- fit_treatment_s
      }
    } else if (treatment_model_type == "joint") {
      treatment_model_args$Y <- A
      treatment_model_args$X <- data.frame(S_dummy, EM_dummy, X)
      fit_treatment <- do.call(what = treatment_model,
                               args = treatment_model_args)
      for (s in unique_S) {
        S_mat <- matrix(0, nrow = n, ncol = no_S - 1)
        colnames(S_mat) <- unique_S[-1]
        if (s %in% unique_S[-1]) S_mat[, s] <- 1
        PrA_XS[, s] <- predict.SuperLearner(fit_treatment,
                                            newdata = data.frame(S_mat, EM_dummy, X))$pred
      }
    } else {
      stop("The 'treatment_model_type' has to be either 'separate' or 'joint'.")
    }

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

    phi_array <- phi_se_array <- array(dim = c(no_S, 3, no_EM),
                                       dimnames = list(unique_S,
                                                       c("A = 1", "A = 0", "Difference"),
                                                       unique_EM))
    for (em in unique_EM) {
      phi <- phi_var <- matrix(nrow = no_S, ncol = 2)
      rownames(phi) <- rownames(phi_var) <- unique_S
      for (s in unique_S) {
        tmp1 <- tmp2 <- matrix(0, nrow = n, ncol = 2)
        I_xs <- which((EM == em) & (S == s))
        kappa <- 1/(length(I_xs)/n)
        tmp1[I_xs, 1] <- pred_Y1[I_xs]
        tmp1[I_xs, 2] <- pred_Y0[I_xs]

        I_xa1 <- which((EM == em) & (A == 1))
        I_xa0 <- which((EM == em) & (A == 0))

        qs <- PrS_X[, s]

        tmp2[I_xa1, 1] <- qs[I_xa1]/eta1[I_xa1] * (Y[I_xa1] - pred_Y1[I_xa1])
        tmp2[I_xa0, 2] <- qs[I_xa0]/eta0[I_xa0] * (Y[I_xa0] - pred_Y0[I_xa0])

        tmp <- tmp1 + tmp2

        phi[s, ] <- kappa * colMeans(tmp)

        tmp1[I_xs, 1] <- pred_Y1[I_xs] - phi[s, 1]
        tmp1[I_xs, 2] <- pred_Y0[I_xs] - phi[s, 2]

        phi_var[s, ] <- kappa/n^2 * colSums((tmp1 + tmp2)^2)
      }
      phi_array[, , em] <- cbind(phi, phi[, 1] - phi[, 2])
      phi_se_array[, , em] <- sqrt(cbind(phi_var, phi_var[, 1] + phi_var[, 2]))
    }
  }

  qt <- qnorm(p = 0.975)
  qtmax <- quantile(apply(abs(matrix(rnorm(no_EM * 1e6),
                                     nrow = no_EM, ncol = 1e6)), 2, max), 0.95)

  lb <- phi_array - qt * phi_se_array
  ub <- phi_array + qt * phi_se_array

  lb_scb <- phi_array - qtmax * phi_se_array
  ub_scb <- phi_array + qtmax * phi_se_array

  # Put results in a data frame
  df_dif <- df_A1 <- df_A0 <-
    data.frame(Source = rep(unique_S, each = no_EM),
               Subgroup = rep(unique_EM, times = no_S),
               Estimate = NA,
               SE = NA,
               ci.lb = NA,
               ci.ub = NA,
               scb.lb = NA,
               scb.ub = NA)

  for (s in unique_S) {
    for (em in unique_EM) {
      rowid <- which(df_dif$Source == s & df_dif$Subgroup == em)
      df_dif$Estimate[rowid] <- phi_array[s, "Difference", em]
      df_dif$SE[rowid] <- phi_se_array[s, "Difference", em]
      df_dif$ci.lb[rowid] <- lb[s, "Difference", em]
      df_dif$ci.ub[rowid] <- ub[s, "Difference", em]
      df_dif$scb.lb[rowid] <- lb_scb[s, "Difference", em]
      df_dif$scb.ub[rowid] <- ub_scb[s, "Difference", em]

      df_A1$Estimate[rowid] <- phi_array[s, "A = 1", em]
      df_A1$SE[rowid] <- phi_se_array[s, "A = 1", em]
      df_A1$ci.lb[rowid] <- lb[s, "A = 1", em]
      df_A1$ci.ub[rowid] <- ub[s, "A = 1", em]
      df_A1$scb.lb[rowid] <- lb_scb[s, "A = 1", em]
      df_A1$scb.ub[rowid] <- ub_scb[s, "A = 1", em]

      df_A0$Estimate[rowid] <- phi_array[s, "A = 0", em]
      df_A0$SE[rowid] <- phi_se_array[s, "A = 0", em]
      df_A0$ci.lb[rowid] <- lb[s, "A = 0", em]
      df_A0$ci.ub[rowid] <- ub[s, "A = 0", em]
      df_A0$scb.lb[rowid] <- lb_scb[s, "A = 0", em]
      df_A0$scb.ub[rowid] <- ub_scb[s, "A = 0", em]
    }
  }

  if (cross_fitting) {
    fit_outcome <- fit_source <- fit_treatment <- NULL
  }

  res <- list(
    df_dif = df_dif,
    df_A0 = df_A0,
    df_A1 = df_A1,
    fit_outcome = fit_outcome,
    fit_source = fit_source,
    fit_treatment = fit_treatment,
    outcome_model_args = outcome_model_args,
    source_model = source_model,
    treatment_model_args = treatment_model_args
  )

  source_names <- unique_S
  subgroup_names <- unique_EM

  res$source_names <- source_names
  res$subgroup_names <- subgroup_names

  class(res) <- 'STE_nested'

  return(res)
}


