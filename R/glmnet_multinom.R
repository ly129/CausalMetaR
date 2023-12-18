#

MN.glmnet <- function(Y, X, newX, s = "lambda.min", family = NULL, obsWeights = NULL, id = NULL, ...) {
  # load required packages
  # require('glmnet')
  Y <- as.factor(Y)
  X <- as.matrix(X)
  if (missing(newX)) {
    newX <- X
  } else {
    newX <- as.matrix(newX)
  }
  fit.MN.glmnet <- glmnet::cv.glmnet(
    x = X,
    y = Y,
    family = "multinomial",
    type.measure = "class",
    weights = obsWeights,
    ...
  )

  # pred is the predicted responses for newX (on the scale of the outcome)
  pred <- predict(
    fit.MN.glmnet,
    newx = newX,
    s = s,
    type = "response"
  )
  # fit returns all objects needed for predict.SL.template
  fit <- list(object = fit.MN.glmnet)
  # return a list with pred and fit
  out <- list(pred = pred[, , 1], fit = fit)
  return(out)
}










#
#
# S_sm <- as.factor(S_sm)
# X_sm <- as.matrix(X_sm)
# if (missing(newX)) {
#   newX <- X_test
# }
# fit.MN.glmnet <- glmnet::cv.glmnet(
#   x = X_sm,
#   y = S_sm,
#   family = "multinomial",
#   type.measure = "class"
# )
#
# # pred is the predicted responses for newX (on the scale of the outcome)
# pred <- predict(
#   fit.MN.glmnet,
#   newx = newX,
#   s = s,
#   type = "response"
# )
# # fit returns all objects needed for predict.SL.template
# fit <- list(object = fit.MN.glmnet)
# # return a list with pred and fit
# out <- list(pred = pred[, , 1], fit = fit)
