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

  return(list(Estimates = phi,
              Variances = phi_var,
              CI_LB = lb,
              CI_UB = ub))#,
              # SCB_LB = lb_scb,
              # SCB_UB = ub_scb))
}


