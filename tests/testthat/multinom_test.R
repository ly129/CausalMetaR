X <- matrix(rnorm(50000), 10000, 5)

# coefficients for each choice
coef1 <- rep(0, 5)
coef2 <- rep(0.5, 5)
coef3 <- rep(-1, 5)

# vector of probabilities
vProb = cbind(exp(X%*%coef1), exp(X%*%coef2), exp(X%*%coef3))

# multinomial draws
mChoices <- t(apply(vProb, 1, rmultinom, n = 1, size = 1))
dfM <- cbind.data.frame(y = apply(mChoices, 1, function(x) which(x == 1)), X)

y <- as.numeric(dfM$y)
x <- dfM[, c(2:6)]
rm(list = "X")

library(SuperLearner)
library(glmnet)
library(nnet)

multinom.fit <- nnet::multinom(factor(y) ~ as.matrix(x) - 1)
# multinom.pred <- predict(multinom.fit, type = "probs")

SL.nnet.multinom <- function(Y, X, newX, family, obsWeights, id, ...) {
  # load required packages
  require('nnet')
  Y <- as.factor(Y)
  if (missing(newX)) {
    newX <- X
  }
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  fit.nnet.multinom <- nnet::multinom(formula = Y ~ X - 1, weights = obsWeights, ...)

  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <- predict(fit.nnet.multinom, newdata = newX, type = "probs")
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = fit.nnet.multinom)
  # declare class of fit for predict.SL.template
  class(fit) <- 'SL.nnet.multinom'
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.nnet.multinom(Y = y, X = x, newX = x, obsWeights = rep(1, 10000))
test <- SuperLearner(Y = y, X = x, newX = x, SL.library = "SL.nnet.multinom", cvControl = list(V = 5))
