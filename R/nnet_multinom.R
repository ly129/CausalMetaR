#

MN.nnet <- function(Y, X, newX, family = NULL, obsWeights = NULL, id = NULL, ...) {
  Y <- as.factor(Y)
  if (missing(newX)) {
    newX <- X
  }

  dtmp <- data.frame(Y = Y, X)

  fit.MN.nnet <- nnet::multinom(formula = Y ~ ., data = dtmp, weights = obsWeights,
                                trace = FALSE, ...)


  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <- predict(fit.MN.nnet, newdata = newX, type = "probs")
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = fit.MN.nnet)
  # return a list with pred and fit
  out <- list(pred = pred, fit = fit)
  return(out)
}
