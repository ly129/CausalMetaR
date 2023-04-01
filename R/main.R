CMetafoR <- function(
    source_model = "S",
    source_model_args,
    treatment_model = "CV.SuperLearner",
    treatment_model_args,
    outcome_model = "CV.SuperLearner",
    outcome_model_args,
) {


}

modelchoices <- function(model1,
                         model2,
                         model2.cv,
                         args1,
                         args2,
                         args2.cv) {
  fit1 <- do.call(model1, args1)
  fit2.cv <- do.call(model2.cv, args2.cv)
  fit2 <- do.call(model2, args2)

  return(list(fit1, fit2))
}
