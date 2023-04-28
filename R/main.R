CMetafoR <- function(
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
    outcome_model_args = list()
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

  return(list(PrS_X,
              PrA_XS,
              predY_AX))
}

