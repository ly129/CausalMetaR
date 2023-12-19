# treatment model type
tmt <- "joint"
# tmt <- "separate"

# source model
srcm <- "MN.nnet"
# srcm <- "MN.glmnet"

# cross_fitting replications
cf.rep <- 5


### STE_nested test
# regular
sn <- STE_nested(
  X = dat_nested[, 2:10],
  Y = dat_nested$Y,
  EM = dat_nested$EM,
  S = dat_nested$S,
  A = dat_nested$A,
  cross_fitting = FALSE,
  source_model = srcm,
  source_model_args = list(trace = FALSE),
  treatment_model_type = tmt,
  treatment_model_args = list(
    family = binomial(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  ),
  outcome_model_args = list(
    family = gaussian(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  )
)
print(sn)
summary(sn)
plot(sn, use_scb = TRUE)
plot(sn, use_scb = FALSE)

# cross-fitting
sncf <- STE_nested(
  X = dat_nested[, 2:10],
  Y = dat_nested$Y,
  EM = dat_nested$EM,
  S = dat_nested$S,
  A = dat_nested$A,
  cross_fitting = TRUE,
  replications = cf.rep,
  source_model = srcm,
  source_model_args = list(trace = FALSE),
  treatment_model_type = tmt,
  treatment_model_args = list(
    family = binomial(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  ),
  outcome_model_args = list(
    family = gaussian(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  )
)
print(sncf)
summary(sncf)
plot(sncf, use_scb = TRUE)
plot(sncf, use_scb = FALSE)


### ATE_nested test
# regular
an <- ATE_nested(
  X = dat_nested[, 2:10],
  Y = dat_nested$Y,
  S = dat_nested$S,
  A = dat_nested$A,
  source_model = srcm,
  source_model_args = list(trace = FALSE),
  treatment_model_type = tmt,
  treatment_model_args = list(
    family = binomial(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  ),
  outcome_model_args = list(
    family = gaussian(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  )
)
summary(an)
print(an)
plot(an)

ancf <- ATE_nested(
  X = dat_nested[, 2:10],
  Y = dat_nested$Y,
  S = dat_nested$S,
  A = dat_nested$A,
  cross_fitting = TRUE,
  replications = cf.rep,
  source_model = srcm,
  source_model_args = list(trace = FALSE),
  treatment_model_type = tmt,
  treatment_model_args = list(
    family = binomial(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  ),
  outcome_model_args = list(
    family = gaussian(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  )
)
summary(ancf)
print(ancf)
plot(ancf)



### STE_external test
se <- STE_external(
  X = dat_nested[, 2:10],
  Y = dat_nested$Y,
  EM = dat_nested$EM,
  S = dat_nested$S,
  A = dat_nested$A,
  X_external = dat_external[, 2:10],
  EM_external = dat_external$EM,
  cross_fitting = FALSE,
  source_model = srcm,
  source_model_args = list(trace = FALSE),
  treatment_model_type = tmt,
  treatment_model_args = list(
    family = binomial(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  ),
  external_model_args = list(
    family = binomial(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  ),
  outcome_model_args = list(
    family = gaussian(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  )
)
print(se)
summary(se)

secf <- STE_external(
  X = dat_nested[, 2:10],
  Y = dat_nested$Y,
  EM = dat_nested$EM,
  S = dat_nested$S,
  A = dat_nested$A,
  X_external = dat_external[, 2:10],
  EM_external = dat_external$EM,
  cross_fitting = FALSE,
  replications = cf.rep,
  source_model = srcm,
  source_model_args = list(trace = FALSE),
  treatment_model_type = tmt,
  treatment_model_args = list(
    family = binomial(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  ),
  external_model_args = list(
    family = binomial(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  ),
  outcome_model_args = list(
    family = gaussian(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  )
)
print(secf)
summary(secf)

### ATE_external test
ae <- ATE_external(
  X = dat_nested[, 2:10],
  Y = dat_nested$Y,
  S = dat_nested$S,
  A = dat_nested$A,
  X_external = dat_external[, 2:10],
  source_model = srcm,
  source_model_args = list(trace = FALSE),
  treatment_model_type = tmt,
  treatment_model_args = list(
    family = binomial(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  ),
  external_model_args = list(
    family = binomial(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  ),
  outcome_model_args = list(
    family = gaussian(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  )
)
summary(ae)
print(ae)

aecf <- ATE_external(
  X = dat_nested[, 2:10],
  Y = dat_nested$Y,
  S = dat_nested$S,
  A = dat_nested$A,
  X_external = dat_external[, 2:10],
  cross_fitting = TRUE,
  replications = cf.rep,
  source_model = srcm,
  source_model_args = list(trace = FALSE),
  treatment_model_type = tmt,
  treatment_model_args = list(
    family = binomial(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  ),
  external_model_args = list(
    family = binomial(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  ),
  outcome_model_args = list(
    family = gaussian(),
    SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
    cvControl = list(V = 5L)
  )
)
summary(aecf)
print(aecf)
