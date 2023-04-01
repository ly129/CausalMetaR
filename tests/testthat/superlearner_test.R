# superlearner

############################
# Setup example dataset.

# Load a dataset from the MASS package.
data(Boston, package = "MASS")

# Extract our outcome variable from the dataframe.
outcome = Boston$medv

# Create a dataframe to contain our explanatory variables.
data = subset(Boston, select = -medv)

# Check structure of our dataframe.
str(data)


# Set a seed for reproducibility in this random sampling.
set.seed(1)

# Reduce to a dataset of 150 observations to speed up model fitting.
train_obs = sample(nrow(data), 150)

# X is our training sample.
x_train = data[train_obs, ]

# Create a holdout set for evaluating model performance.
# Note: cross-validation is even better than a single holdout sample.
x_holdout = data[-train_obs, ]

# Create a binary outcome variable: towns in which median home value is > 22,000.
outcome_bin = as.numeric(outcome > 22)

y_train = outcome_bin[train_obs]
y_holdout = outcome_bin[-train_obs]

# Review the outcome variable distribution.
table(y_train, useNA = "ifany")


library(SuperLearner)

listWrappers()


glm(y_train ~ as.matrix(x_train), family = binomial())

sl.binom.glm = SuperLearner(Y = y_train,
                            X = x_train,
                            family = binomial(),
                            SL.library = "SL.glm")
sl.binom.glm$fitLibrary$SL.glm_All$object$coefficients


glmnet.5.fold <- create.Learner("SL.glmnet", params = list(nfolds = 5))
sl.glmnet.5.fold <- glmnet.5.fold$names
sl.binom.l1 <- SuperLearner(Y = y_train,
                            X = x_train,
                            family = binomial(),
                            SL.library = "SL.glmnet",
                            cvControl = list(V = 5))
predict(sl.binom.l1, newdata = x_train)


cv_sl = CV.SuperLearner(Y = y_train, X = x_train, family = binomial(),
                        # For a real analysis we would use V = 10.
                        cvControl = list(V = 5), innerCvControl = list(list(V=10)),
                        SL.library = c("SL.mean", "SL.glmnet", "SL.glm", "SL.xgboost"))
plot(cv_sl)
predict(cv_sl$AllSL)

