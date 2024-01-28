#' Estimating the Average Treatment Effect (ATE) in an external target population using multi-source data
#'
#' @description
#' Doubly-robust and efficient estimator for the ATE in an external target population using multi-source data.
#'
#' @param X Data frame (or matrix) containing the covariate data in the multi-source data. It should have \eqn{n} rows and \eqn{p} columns. Character variables will be converted to factors.
#' @param X_external Data frame (or matrix) containing the covariate data in the external target population. It should have \eqn{n_0} rows and \eqn{p} columns. This is the external data counterpart to the \code{X} argument.
#' @param Y Vector of length \eqn{n} containing the outcome.
#' @param S Vector of length \eqn{n} containing the source indicator. If \code{S} is a factor, it will maintain its level order; otherwise it will be converted to a factor with the default level order. The order will be carried over to the outputs and plots.
#' @param A Vector of length \eqn{n} containing the binary treatment (1 for treated and 0 for untreated).
#' @param cross_fitting Logical specifying whether sample splitting and cross fitting should be used.
#' @param replications Integer specifying the number of sample splitting and cross fitting replications to perform, if \code{cross_fitting = TRUE}. The default is \code{10L}.
#' @param source_model Character string specifying the (penalized) multinomial logistic regression for estimating the source model. It has two options: "\code{MN.glmnet}" (default) and "\code{MN.nnet}", which use \pkg{glmnet} and \pkg{nnet} respectively.
#' @param source_model_args List specifying the arguments for the source model (in \pkg{glmnet} or \pkg{nnet}).
#' @param treatment_model_type Character string specifying how the treatment model is estimated. Options include "\code{separate}" (default) and "\code{joint}". If "\code{separate}", the treatment model (i.e., \eqn{P(A=1|X, S=s)}) is estimated by regressing \eqn{A} on \eqn{X} within each specific internal population \eqn{S=s}. If "\code{joint}", the treatment model is estimated by regressing \eqn{A} on \eqn{X} and \eqn{S} using the multi-source population.
#' @param treatment_model_args List specifying the arguments for the treatment model (in \pkg{SuperLearner}).
#' @param external_model_args List specifying the arguments for the external model (in \pkg{SuperLearner}).
#' @param outcome_model_args List specifying the arguments for the outcome model  (in \pkg{SuperLearner}).
#' @param show_progress Logical specifying whether to print a progress bar for the cross-fit replicates completed, if \code{cross_fitting = TRUE}.
#'
#' @details
#'
#' \strong{Data structure:}
#'
#' The multi-source dataset consists the outcome \code{Y}, source \code{S}, treatment \code{A}, and covariates \code{X} (\eqn{n \times p}) in the internal populations. The data sources can be trials, observational studies, or a combination of both.
#'
#' The external dataset contains only covariates \code{X_external} (\eqn{n_0 \times p}).
#'
#' \strong{Estimation of nuissance parameters:}
#'
#' The following models are fit:
#' \itemize{
#' \item External model: \eqn{q(X)=P(R=1|X)}, where \eqn{R} takes value 1 if the subject belongs to any of the internal dataset and 0 if the subject belongs to the external dataset
#' \item Propensity score model: \eqn{\eta_a(X)=P(A=a|X)}. We perform the decomposition \eqn{P(A=a|X)=\sum_{s} P(A=a|X, S=s)P(S=s|X)} and estimate \eqn{P(A=1|X, S=s)} (i.e., the treatment model) and \eqn{P(S=s|X)} (i.e., the source model).
#' \item Outcome model: \eqn{\mu_a(X)=E(Y|X, A=a)}
#' }
#' The models are estimated by \pkg{SuperLearner} with the exception of the source model which is estimated by \pkg{glmnet} or \pkg{nnet}.
#'
#' \strong{ATE estimation:}
#'
#' The ATE estimator is
#' \deqn{
#'  \dfrac{\widehat \kappa}{N}\sum\limits_{i=1}^{N} \Bigg[ I(R_i = 0) \widehat \mu_a(X_i)
#'  +I(A_i = a, R_i=1) \dfrac{1-\widehat q(X_i)}{\widehat \eta_a(X_i)\widehat q(X_i)}  \Big\{ Y_i - \widehat \mu_a(X_i) \Big\} \Bigg],
#' }
#' where \eqn{N=n+n_0}, and \eqn{\widehat \kappa=\{N^{-1} \sum_{i=1}^N I(R_i=0)\}^{-1}}.
#'
#' This estimator is doubly robust and non-parametrically efficient.
#' To achieve non-parametric efficiency and asymptotic normality, it requires that \eqn{||\widehat \mu_a(X) -\mu_a(X)||\big\{||\widehat \eta_a(X) -\eta_a(X)||+||\widehat q(X) -q(X)||\big\}=o_p(n^{-1/2})}.
#' In addition, sample splitting and cross-fitting can be performed to avoid the Donsker class assumption.
#'
#' When a data source is a randomized trial, it is still recommended to estimate the propensity score for optimal efficiency.
#'
#' @return An object of class "ATE_external". This object is a list with the following elements:
#'   \item{df_dif}{A data frame containing the treatment effect (mean difference) estimates for the extenal data.}
#'   \item{df_A0}{A data frame containing the potential outcome mean estimates under A = 0 for the extenal data.}
#'   \item{df_A1}{A data frame containing the potential outcome mean estimates under A = 1 for the extenal data.}
#'   \item{fit_outcome}{Fitted outcome model.}
#'   \item{fit_source}{Fitted source model.}
#'   \item{fit_treatment}{Fitted treatment model(s).}
#'   \item{fit_external}{Fitted external model.}
#'
#' @references Dahabreh, I.J., Robertson, S.E., Petito, L.C., Hernán, M.A. and Steingrimsson, J.A.. (2019) \emph{Efficient and robust methods for causally
#' interpretable meta‐analysis: Transporting inferences from multiple randomized trials to a target population}, Biometrics.
#'
#' @examples
#' \donttest{
#' ae <- ATE_external(
#'   X = dat_multisource[, 1:10],
#'   Y = dat_multisource$Y,
#'   S = dat_multisource$S,
#'   A = dat_multisource$A,
#'   X_external = dat_external[, 1:10],
#'   source_model = "MN.glmnet",
#'   source_model_args = list(),
#'   treatment_model_type = "separate",
#'   treatment_model_args = list(
#'     family = binomial(),
#'     SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
#'     cvControl = list(V = 5L)
#'   ),
#'   external_model_args = list(
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
#' }
#'
#' @import metafor
#' @import SuperLearner
#' @importFrom stats model.matrix predict qnorm quantile rnorm median
#' @importFrom graphics text
#'
#' @export

ATE_external <- function(
    X,
    X_external,
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
    outcome_model_args = list(),
    show_progress = TRUE
) {
  # For future possibilities
  treatment_model <- external_model <- outcome_model <- "SuperLearner"

  # Total sample size
  n1 <- nrow(X)
  n0 <- nrow(X_external)
  n <- n0 + n1

  # Checking lengths of all data inputs
  if (length(Y) != n1 | length(S) != n1 | length(A) != n1){
    stop(paste0('All data inputs for the sources must have the same length:\n',
                'X has length ', n1, '.\n',
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

  # Covariates
  if (ncol(X) != ncol(X_external)){
    stop(paste0('X and X_external must have the same number of columns.'))
  }
  # Converting character variables into factor variables
  X <- as.data.frame(sapply(X, FUN = function(xx) {
    if (is.character(xx)) {
      factor(xx)
    } else {
      xx
    }
  }))
  X_external <- as.data.frame(sapply(X_external, FUN = function(xx) {
    if (is.character(xx)) {
      factor(xx)
    } else {
      xx
    }
  }))

  # Converting factor variables into dummy variables
  X <- data.frame(model.matrix(~ ., data = X)[, -1])
  X_external <- data.frame(model.matrix(~ ., data = X_external)[, -1])

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
    K <- 5L
    psi_array_cf <- psi_se_array_cf <- array(dim = c(3, K, replications),
                                             dimnames = list(c("A = 1", "A = 0", "Difference"),
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
        S_test <- S[test_id]

        X_sm <- X[sm_id, ]
        Y_sm <- Y[sm_id]
        A_sm <- A[sm_id]
        S_sm <- S[sm_id]

        X_tm <- X[tm_id, ]
        Y_tm <- Y[tm_id]
        A_tm <- A[tm_id]
        S_tm <- S[tm_id]
        if (treatment_model_type == "joint") S_dummy_tm <- S_dummy[tm_id, ]

        X_om <- X[om_id, ]
        Y_om <- Y[om_id]
        A_om <- A[om_id]
        S_om <- S[om_id]

        X_xm <- X[xm_id, ]
        Y_xm <- Y[xm_id]
        A_xm <- A[xm_id]
        S_xm <- S[xm_id]

        X_ext_xm <- X_external[ext_xm_id, ]

        # source model
        if (source_model %in% c("MN.glmnet", "MN.nnet")) {
          source_model_args$Y <- S
          source_model_args$X <- X
          source_model_args$newX <- X_test
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

        # external model
        external_model_args$Y <- c(rep(1, length(xm_id)), rep(0, length(ext_xm_id)))
        external_model_args$X <- rbind(X_xm, X_ext_xm)
        fit_external <- do.call(what = external_model, args = external_model_args)
        PrR_X <- predict.SuperLearner(fit_external, newdata = X_test)$pred

        # outcome model
        outcome_model_args$Y <- Y_om
        outcome_model_args$X <- data.frame(A = A_om, X_om)
        fit_outcome <- do.call(what = outcome_model, args = outcome_model_args)
        pred_Y1 <- predict.SuperLearner(fit_outcome,
                                        newdata = data.frame(A = 1, X_test))$pred
        pred_Y0 <- predict.SuperLearner(fit_outcome,
                                        newdata = data.frame(A = 0, X_test))$pred
        predY_AX <- cbind(pred_Y1, pred_Y0)

        # estimators

        eta1 <- rowSums(PrA_XS * PrS_X)
        eta0 <- rowSums((1 - PrA_XS) * PrS_X)

        pred_Y1_X0 <- c(predict.SuperLearner(fit_outcome,
                                             newdata = data.frame(A = 1, X_ext_xm))$pred)
        pred_Y0_X0 <- c(predict.SuperLearner(fit_outcome,
                                             newdata = data.frame(A = 0, X_ext_xm))$pred)
        pred_Y1_X1 <- c(predict.SuperLearner(fit_outcome,
                                             newdata = data.frame(A = 1, X_test))$pred)
        pred_Y0_X1 <- c(predict.SuperLearner(fit_outcome,
                                             newdata = data.frame(A = 0, X_test))$pred)

        gamma <- n/n0 # length(I_xr)

        tmp1 <- matrix(0, nrow = length(ext_xm_id), ncol = 2)
        tmp1[, 1] <- pred_Y1_X0   #[I_xr, ]
        tmp1[, 2] <- pred_Y0_X0   #[I_xr, ]

        tmp2 <- matrix(0, nrow = length(test_id), ncol = 2)
        I_xa1 <- which(A_test == 1)
        I_xa0 <- which(A_test == 0)

        tmp2[I_xa1, 1] <- (1 - PrR_X[I_xa1])/PrR_X[I_xa1]/eta1[I_xa1] * (Y_test[I_xa1] - pred_Y1_X1[I_xa1])
        tmp2[I_xa0, 2] <- (1 - PrR_X[I_xa0])/PrR_X[I_xa0]/eta0[I_xa0] * (Y_test[I_xa0] - pred_Y0_X1[I_xa0])

        tmp <- colSums(rbind(tmp1, tmp2))
        psi <- gamma/(length(ext_xm_id)+length(test_id)) * tmp

        tmp1 <- tmp1 - rep(psi, each = length(ext_xm_id))
        psi_var <- gamma^2/(length(ext_xm_id)+length(test_id))^2 * colSums(rbind(tmp1, tmp2)^2)

        psi_array_cf[, k, r] <- c(psi, psi[1] - psi[2])
        psi_se_array_cf[, k, r] <- sqrt(c(psi_var, psi_var[1] + psi_var[2]))
      } # end of k loop
      if (show_progress){
        pb$tick()
      }
    } # end of r loop

    psi_array <- apply(apply(psi_array_cf, MARGIN = c(1, 3), FUN = mean), MARGIN = 1, FUN = median)
    psi_se_array <- apply(apply(psi_se_array_cf, MARGIN = c(1, 3), FUN = mean), MARGIN = 1, FUN = median)
    # end of ATE_external with cross-fitting
  } else {
    # start of regular STE_external

    if (source_model %in% c("MN.glmnet", "MN.nnet")) {
      source_model_args$Y <- S
      source_model_args$X <- X
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
        S_mat <- matrix(0, nrow = n1, ncol = no_S - 1)
        colnames(S_mat) <- unique_S[-1]
        if (s %in% unique_S[-1]) S_mat[, s] <- 1
        PrA_XS[, s] <- predict.SuperLearner(fit_treatment,
                                            newdata = data.frame(S_mat, X))$pred
      }
    } else {
      stop("The 'treatment_model_type' has to be either 'separate' or 'joint'.")
    }

    external_model_args$Y <- c(rep(1, n1), rep(0, n0))
    external_model_args$X <- rbind(X, X_external)
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
    eta1 <- rowSums(PrA_XS * PrS_X)
    eta0 <- rowSums((1 - PrA_XS) * PrS_X)

    pred_Y1_X0 <- predict.SuperLearner(fit_outcome,
                                       newdata = data.frame(A = 1, X_external))$pred
    pred_Y0_X0 <- predict.SuperLearner(fit_outcome,
                                       newdata = data.frame(A = 0, X_external))$pred
    pred_Y1_X1 <- predict.SuperLearner(fit_outcome,
                                       newdata = data.frame(A = 1, X))$pred
    pred_Y0_X1 <- predict.SuperLearner(fit_outcome,
                                       newdata = data.frame(A = 0, X))$pred

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
    psi <- gamma/n * tmp

    tmp1 <- tmp1 - rep(psi, each = n0)
    psi_var <- gamma^2/n^2 * colSums(rbind(tmp1, tmp2)^2)

    psi_array <- c(psi, psi[1] - psi[2])
    psi_se_array <- sqrt(c(psi_var, psi_var[1] + psi_var[2]))

    names(psi_array) <- names(psi_se_array) <- c("A = 1", "A = 0", "Difference")
  }

  qt <- qnorm(p = 0.975)

  lb <- psi_array - qt * sqrt(psi_se_array)
  ub <- psi_array + qt * sqrt(psi_se_array)

  df_dif <-
    data.frame(Estimate = psi_array["Difference"],
               SE = psi_se_array["Difference"],
               ci.lb = lb["Difference"],
               ci.ub = ub["Difference"])
  df_A0 <-
    data.frame(Estimate = psi_array["A = 0"],
               SE = psi_se_array["A = 0"],
               ci.lb = lb["A = 0"],
               ci.ub = ub["A = 0"])
  df_A1 <-
    data.frame(Estimate = psi_array["A = 1"],
               SE = psi_se_array["A = 1"],
               ci.lb = lb["A = 1"],
               ci.ub = ub["A = 1"])

  if (cross_fitting) {
    fit_outcome <- fit_source <- fit_treatment <- fit_external <- NULL
  }

  res <- list(df_dif = df_dif,
              df_A0 = df_A0,
              df_A1 = df_A1,
              fit_outcome = fit_outcome,
              fit_source = fit_source,
              fit_treatment = fit_treatment,
              fit_external = fit_external,
              outcome_model_args = outcome_model_args,
              source_model = source_model,
              treatment_model_args = treatment_model_args,
              external_model_args = external_model_args)

  class(res) <- 'ATE_external'

  return(res)
}


