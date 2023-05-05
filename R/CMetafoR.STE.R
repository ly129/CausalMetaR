CMetafoR.STE <- function(
    X,
    Y,
    S, # integer sequence starting from 1
    A,
    source_model = "SL.glmnet.multinom",
    source_model_args = list(),
    treatment_model_type = "separate",
    treatment_model = "SuperLearner",
    treatment_model_args = list(),
    outcome_model = "SuperLearner",
    outcome_model_args = list(),
    x_tilde
) {
  # Total sample size
  n <- nrow(X)

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

  PrA_XS <- matrix(nrow = n, ncol = no_S)
  if (treatment_model_type == "separate") {
    for (s in 1:no_S) {
      id_s <- which(S == s)
      treatment_model_args$Y <- A[id_s]
      treatment_model_args$X <- X[id_s, ]
      fit_treatment_s <- do.call(what = treatment_model,
                                 args = treatment_model_args)
      PrA_XS[, s] <- predict.SuperLearner(fit_treatment_s, newdata = X)$pred
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
    I_xs <- which((X[, 1] == x_tilde) & (S == s))
    kappa <- 1/(length(I_xs)/n)
    tmp1[I_xs, 1] <- pred_Y1[I_xs]
    tmp1[I_xs, 2] <- pred_Y0[I_xs]

    I_xa1 <- which((X[, 1] == x_tilde) & (A == 1))
    I_xa0 <- which((X[, 1] == x_tilde) & (A == 0))

    qs <- PrS_X[, s]

    tmp2[I_xa1, 1] <- qs[I_xa1]/eta1[I_xa1] * (Y[I_xa1] - pred_Y1[I_xa1])
    tmp2[I_xa0, 2] <- qs[I_xa0]/eta0[I_xa0] * (Y[I_xa0] - pred_Y0[I_xa0])

    tmp <- tmp1 + tmp2

    psi[s, ] <- kappa * colMeans(tmp)

    tmp1[I_xs, 1] <- pred_Y1[I_xs] - psi[s, 1]
    tmp1[I_xs, 2] <- pred_Y0[I_xs] - psi[s, 2]

    psi_var[s, ] <- kappa/n^2 * colSums((tmp1 + tmp2)^2)
  }

  rownames(psi) <- paste0("S=", unique_S)
  colnames(psi) <- paste0("A=", c(1, 0))

  rownames(psi_var) <- paste0("S=", unique_S)
  colnames(psi_var) <- paste0("A=", c(1, 0))

  lb <- psi - qnorm(p = 0.975) * sqrt(psi_var)
  ub <- psi + qnorm(p = 0.975) * sqrt(psi_var)

  tmax <- apply(abs(matrix(rnorm(length(unique(X[,1])) * 1e6),
                           nrow = length(unique(X[,1])), ncol = 1e6)), 2, max)
  qtmax <- quantile(tmax, 0.95)

  lb_scb <- psi - qtmax * sqrt(psi_var)
  ub_scb <- psi + qtmax * sqrt(psi_var)

  return(list(Estimates = psi,
              Variances = psi_var,
              CI_LB = lb,
              CI_UB = ub,
              SCB_LB = lb_scb,
              SCB_UB = ub_scb))
}


