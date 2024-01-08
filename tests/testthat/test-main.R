outcome_model_args <- list(
  family = gaussian(),
  SL.library = c("SL.glmnet", "SL.glm"),
  cvControl = list(V = 5L))

treatment_model_args <- list(
  family = binomial(),
  SL.library = c("SL.glmnet", "SL.glm"),
  cvControl = list(V = 5L))

external_model_args = list(
  family = binomial(),
  SL.library = c("SL.glmnet", "SL.glm"),
  cvControl = list(V = 5L))

replications <- 2

set.seed(1234)
for (cross_fitting in c(FALSE)){
  for (treatment_model_type in c('separate', 'joint')){
    for (source_model in c('MN.nnet', 'MN.glmnet')){

      test_that(
        paste('ATE_internal with: cross_fitting = ', cross_fitting, ', treatment_model_type = ', treatment_model_type, ', source_model = ', source_model), {
          expect_no_error({
              res <- ATE_internal(
                X = dat_multisource[, 1:5],
                Y = dat_multisource$Y,
                S = dat_multisource$S,
                A = dat_multisource$A,
                cross_fitting = cross_fitting,
                source_model = source_model,
                treatment_model_type = treatment_model_type,
                treatment_model_args = treatment_model_args,
                outcome_model_args = outcome_model_args,
                replications = replications)
              capture.output(print(res))
              capture.output(summary(res))
            })
        })

      test_that(
        paste('STE_internal with: cross_fitting = ', cross_fitting, ', treatment_model_type = ', treatment_model_type, ', source_model = ', source_model), {
          expect_no_error({
              res <- STE_internal(
                X = dat_multisource[, 1:5],
                Y = dat_multisource$Y,
                EM = dat_multisource$EM,
                S = dat_multisource$S,
                A = dat_multisource$A,
                cross_fitting = cross_fitting,
                source_model = source_model,
                treatment_model_type = treatment_model_type,
                treatment_model_args = treatment_model_args,
                outcome_model_args = outcome_model_args,
                replications = replications)
              capture.output(print(res))
              capture.output(summary(res))
            })
        })

      test_that(
        paste('ATE_external with: cross_fitting = ', cross_fitting, ', treatment_model_type = ', treatment_model_type, ', source_model = ', source_model), {
          expect_no_error({
              res <- ATE_external(
                X = dat_multisource[, 2:5],
                X_external = dat_external[, 2:5],
                Y = dat_multisource$Y,
                S = dat_multisource$S,
                A = dat_multisource$A,
                cross_fitting = cross_fitting,
                source_model = source_model,
                treatment_model_type = treatment_model_type,
                treatment_model_args = treatment_model_args,
                outcome_model_args = outcome_model_args,
                external_model_args = external_model_args,
                replications = replications)
              capture.output(print(res))
              capture.output(summary(res))
            })
        })

      test_that(
        paste('STE_external with: cross_fitting = ', cross_fitting, ', treatment_model_type = ', treatment_model_type, ', source_model = ', source_model), {
          expect_no_error({
              res <- STE_external(
                X = dat_multisource[, 2:5],
                X_external = dat_external[, 2:5],
                Y = dat_multisource$Y,
                EM = dat_multisource$EM,
                EM_external = dat_external$EM,
                S = dat_multisource$S,
                A = dat_multisource$A,
                cross_fitting = cross_fitting,
                source_model = source_model,
                treatment_model_type = treatment_model_type,
                treatment_model_args = treatment_model_args,
                outcome_model_args = outcome_model_args,
                external_model_args = external_model_args,
                replications = replications)
              capture.output(print(res))
              capture.output(summary(res))
            })
        })
    }
  }
}
