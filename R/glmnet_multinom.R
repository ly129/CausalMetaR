#

glmnet.multinom <- function(Y, X, newX, s = "lambda.min", family = NULL, obsWeights = NULL, id = NULL, ...) {
  # load required packages
  # require('glmnet')
  Y <- as.factor(Y)
  X <- as.matrix(X)
  if (missing(newX)) {
    newX <- X
  }
  fit.glmnet.multinom <- glmnet::cv.glmnet(
    x = X,
    y = Y,
    family = "multinomial",
    type.measure = "class",
    weights = obsWeights,
    ...
  )

  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <- predict(
    fit.glmnet.multinom,
    newx = newX,
    s = s,
    type = "response"
  )
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = fit.glmnet.multinom)
  # return a list with pred and fit
  out <- list(pred = pred[, , 1], fit = fit)
  return(out)
}
