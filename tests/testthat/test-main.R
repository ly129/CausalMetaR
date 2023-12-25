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
        paste('ATE_nested with: cross_fitting = ', cross_fitting, ', treatment_model_type = ', treatment_model_type, ', source_model = ', source_model), {
          expect_no_error({
              res <- ATE_nested(
                X = dat_nested[, 1:5],
                Y = dat_nested$Y,
                S = dat_nested$S,
                A = dat_nested$A,
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
        paste('STE_nested with: cross_fitting = ', cross_fitting, ', treatment_model_type = ', treatment_model_type, ', source_model = ', source_model), {
          expect_no_error({
              res <- STE_nested(
                X = dat_nested[, 1:5],
                Y = dat_nested$Y,
                EM = dat_nested$EM,
                S = dat_nested$S,
                A = dat_nested$A,
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
                X = dat_nested[, 2:5],
                X_external = dat_external[, 2:5],
                Y = dat_nested$Y,
                S = dat_nested$S,
                A = dat_nested$A,
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
                X = dat_nested[, 2:5],
                X_external = dat_external[, 2:5],
                Y = dat_nested$Y,
                EM = dat_nested$EM,
                EM_external = dat_external$EM,
                S = dat_nested$S,
                A = dat_nested$A,
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
