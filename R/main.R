CMetafoR <- function(
    X,
    Y,
    S, # integer sequence starting from 1
    A,
    source_model = "SL.glmnet.multinom",
    source_model_args = list(),
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
    fit_source <- do.call(what = source_model, args = source_model_args)
    pred_S <- fit_source$pred
  } else {
    stop("Currently only support `SL.glmnet.multinom` and `SL.nnet.multinom`.")
  }

  fit_treatment <- do.call(what = treatment_model, args = treatment_model_args)
  pred_A <- predict.SuperLearner(fit_treatment)$pred

  fit_outcome <- do.call(what = outcome_model, args = outcome_model_args)
  pred_Y <- predict.SuperLearner(fit_outcome)

  return(list(pred_S,
              pred_A,
              pred_Y))
}

