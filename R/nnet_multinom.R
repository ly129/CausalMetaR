#

SL.nnet.multinom <- function(Y, X, newX, family = NULL, obsWeights = NULL, id = NULL, ...) {
  # load required packages
  # require('nnet')
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
