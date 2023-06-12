#' Transport ATE from multi-source data to an internal source-specific data
#'
#' @description
#' Estimate the average treatment effect of an internal source-specific target population. The nuisance parameters are allowed to be estimated by SuperLearner.   
#' 
#' @param x Predictor matrix with dimension \eqn{nm * p}, where \eqn{n} is the number of subjects, \eqn{m} is the maximum observation time, and \eqn{p} is the number of predictors. See Details.
#' @param ID The ID of each subjects, each subject has one ID (many rows in \code{x} share one \code{ID}).
#' @param time Represents the start of each time interval.
#' @param time2 Represents the stop of each time interval.
#' @param event Indicator of event. \code{event = 1} when event occurs and \code{event = 0} otherwise.
#' @param lambda Sequence of regularization coefficients \eqn{\lambda}'s.
#' @param group \eqn{G * G} matrix describing the relationship between the groups of variables, where \eqn{G} represents the number of groups. Denote the \eqn{i}-th group of variables by \eqn{g_i}. The \eqn{(i,j)} entry is \code{1} if and only if \eqn{i\neq j} and \eqn{g_i} is a child group (subset) of \eqn{g_j}, and is \code{0} otherwise. See Examples and Details.
#' @param group_variable \eqn{p * G} matrix describing the relationship between the groups and the variables. The \eqn{(i,j)} entry is \code{1} if and only if variable \eqn{i} is in group \eqn{g_j}, but not in any child group of \eqn{g_j}, and is \code{0} otherwise. See Examples and Details.
#' @param penalty_weights Optional, vector of length \eqn{G} specifying the group-specific penalty weights. If not specified, the default value is \eqn{\mathbf{1}_G}. Modify with caution.
#' @param par_init Optional, vector of initial values of the optimization algorithm. Default initial value is zero for all \eqn{p} variables.
#' @param stepsize_init Initial value of the stepsize of the optimization algorithm. Default is 1.
#' @param stepsize_shrink Factor in \eqn{(0,1)} by which the stepsize shrinks in the backtracking linesearch. Default is 0.8.
#' @param tol Convergence criterion. Algorithm stops when the \eqn{l_2} norm of the difference between two consecutive updates is smaller than \code{tol}.
#' @param maxit Maximum number of iterations allowed.
#' @param verbose Logical, whether progress is printed.
#' 
#' @examples 
#' # g3 in g1 -> grp_31 = 1
#' # g3 in g2 -> grp_32 = 1
#' # g5 in g2 -> grp_52 = 1
#' # g5 in g4 -> grp_54 = 1
#' grp <- matrix(c(0, 0, 0, 0, 0,
#'                 0, 0, 0, 0, 0,
#'                 1, 1, 0, 0, 0,
#'                 0, 0, 0, 0, 0,
#'                 0, 1, 0, 1, 0),
#'               ncol = 5, byrow = TRUE)
#'
#' # Variable A1 is in g1 only: grp.var_11 = 1
#' # Variable A1B is in g1 and g3, but g3 is a child group of g1,
#' # so grp.var_63 = 1 while grp.var_61 = 0.
#' grp.var <- matrix(c(1, 0, 0, 0, 0, #A1
#'                     1, 0, 0, 0, 0, #A2
#'                     0, 0, 0, 1, 0, #C1
#'                     0, 0, 0, 1, 0, #C2
#'                     0, 1, 0, 0, 0, #B
#'                     0, 0, 1, 0, 0, #A1B
#'                     0, 0, 1, 0, 0, #A2B
#'                     0, 0, 0, 0, 1, #C1B
#'                     0, 0, 0, 0, 1  #C2B
#'                    ), ncol = 5, byrow = TRUE)
#' eta_g <- rep(1, 5)
#' x <- as.matrix(sim[, c("A1","A2","C1","C2","B",
#'                        "A1B","A2B","C1B","C2B")])
#' lam.seq <- 10^seq(0, -2, by = -0.2)
#' 
#' fit <- sox(x = x,
#'               ID = sim$Id,
#'               time = sim$Start,
#'               time2 = sim$Stop,
#'               event = sim$Event,
#'               lambda = lam.seq,
#'               group = grp,
#'               group_variable = grp.var,
#'               penalty_weights = eta_g,
#'               tol = 1e-4,
#'               maxit = 1e3,
#'               verbose = FALSE)
#' @details
#' The predictor matrix should be of dimension \eqn{nm * p}. Each row records the values of covariates for one subject at one time, for example, the values at the day from \code{time} (Start) to \code{time2} (Stop). An example dataset \code{\link{sim}} is provided. The dataset has the same format produced by the \code{R} package \pkg{PermAlgo}. 
#' The specification of arguments \code{group} and \code{group_variable} for the grouping structure can be found in \url{http://thoth.inrialpes.fr/people/mairal/spams/doc-R/html/doc_spams006.html#sec27}, the same as the grouping structure specification in the \code{R} package \pkg{spams}.
#'
#' In the Examples below, \eqn{p=9,G=5}, the group structure is: \deqn{g_1 = \{A_{1}, A_{2}, A_{1}B, A_{2}B\},} \deqn{g_2  = \{B, A_{1}B, A_{2}B, C_{1}B, C_{2}B\},} \deqn{g_3  = \{A_{1}B, A_{2}B\},} \deqn{g_4  = \{C_1, C_2, C_{1}B, C_{2}B\},} \deqn{g_5  = \{C_{1}B, C_{2}B\}.}
#' 
#' where \eqn{g_3} is a subset of \eqn{g_1} and \eqn{g_2}, and \eqn{g_5} is a subset of \eqn{g_2} and \eqn{g_4}.
#' @return A list with the following three elements.
#'   \item{lambdas}{The user-specified regularization coefficients \code{lambda} sorted in decreasing order.}
#'   \item{estimates}{A matrix, with each column corresponding to the coefficient estimates at each \eqn{\lambda} in \code{lambdas}.}
#'   \item{iterations}{A vector of number of iterations it takes to converge at each \eqn{\lambda} in \code{lambdas}.}

CMetafoR.ATE.S <- function(
    X,
    Y,
    S, # integer sequence starting from 1
    A,
    source_model = "SL.glmnet.multinom",
    source_model_args = list(),
    treatment_model_type = "separate",
    treatment_model = "SuperLearner",
    treatment_model_args = list(),
    outcome_model = "SuperLearner",
    outcome_model_args = list()
) {
  # Total sample size
  n <- nrow(X)

  # Number of sources - with format check
  unique_S <- sort(unique(S))
  no_S <- length(unique_S)
  if (! identical(seq(no_S), unique_S)) {
    stop(paste("Source", setdiff(seq(max(S)), unique_S), "is missing. "))
  }

  if (source_model %in% c("SL.glmnet.multinom", "SL.nnet.multinom")) {
    source_model_args$Y <- S
    source_model_args$X <- X
    fit_source <- do.call(what = source_model, args = source_model_args)
    PrS_X <- fit_source$pred
  } else {
    stop("Currently only support `SL.glmnet.multinom` and `SL.nnet.multinom`.")
  }

  PrA_XS <- matrix(nrow = n, ncol = no_S)
  if (treatment_model_type == "separate") {
    for (s in 1:no_S) {
      id_s <- which(S == s)
      treatment_model_args$Y <- A[id_s]
      treatment_model_args$X <- X[id_s, ]
      fit_treatment_s <- do.call(what = treatment_model,
                                 args = treatment_model_args)
      PrA_XS[, s] <- predict.SuperLearner(fit_treatment_s, newdata = X)$pred
    }
  } else if (treatment_model_type == "joint") {
    treatment_model_args$Y <- A
    treatment_model_args$X <- cbind(X, S)
    fit_treatment <- do.call(what = treatment_model,
                             args = treatment_model_args)
    for (s in 1:no_S) {
      PrA_XS[, s] <- predict.SuperLearner(fit_treatment,
                                          newdata = cbind(X, s))$pred
    }
  } else {
    stop("Type has to be either 'separate' or 'joint'.")
  }

  outcome_model_args$Y <- Y
  outcome_model_args$X <- data.frame(A, X)
  fit_outcome <- do.call(what = outcome_model, args = outcome_model_args)
  pred_Y1 <- predict.SuperLearner(fit_outcome,
                                  newdata = data.frame(A = 1, X))$pred
  pred_Y0 <- predict.SuperLearner(fit_outcome,
                                  newdata = data.frame(A = 0, X))$pred
  predY_AX <- cbind(pred_Y1, pred_Y0)

  # estimators

  eta1 <- PrA_XS * PrS_X
  eta0 <- (1 - PrA_XS) * PrS_X

  psi <- psi_var <- matrix(nrow = no_S, ncol = 2)

  for (s in unique_S) {
    tmp1 <- tmp2 <- matrix(0, nrow = n, ncol = 2)
    I_xs <- which(S == s)
    kappa <- 1/(length(I_xs)/n)
    tmp1[I_xs, 1] <- pred_Y1[I_xs]
    tmp1[I_xs, 2] <- pred_Y0[I_xs]

    I_xa1 <- which(A == 1)
    I_xa0 <- which(A == 0)

    qs <- PrS_X[, s]

    tmp2[I_xa1, 1] <- qs[I_xa1]/eta1[I_xa1] * (Y[I_xa1] - pred_Y1[I_xa1])
    tmp2[I_xa0, 2] <- qs[I_xa0]/eta0[I_xa0] * (Y[I_xa0] - pred_Y0[I_xa0])

    tmp <- tmp1 + tmp2

    psi[s, ] <- kappa * colMeans(tmp)

    tmp1[I_xs, 1] <- pred_Y1[I_xs] - psi[s, 1]
    tmp1[I_xs, 2] <- pred_Y0[I_xs] - psi[s, 2]

    psi_var[s, ] <- kappa/n^2 * colSums((tmp1 + tmp2)^2)
  }

  rownames(psi) <- paste0("S=", unique_S)
  colnames(psi) <- paste0("A=", c(1, 0))

  rownames(psi_var) <- paste0("S=", unique_S)
  colnames(psi_var) <- paste0("A=", c(1, 0))

  lb <- psi - qnorm(p = 0.975) * sqrt(psi_var)
  ub <- psi + qnorm(p = 0.975) * sqrt(psi_var)

  # tmax <- apply(abs(matrix(rnorm(length(unique(X[,1])) * 1e6),
  #                          nrow = length(unique(X[,1])), ncol = 1e6)), 2, max)
  # qtmax <- quantile(tmax, 0.95)
  #
  # lb_scb <- psi - qtmax * sqrt(psi_var)
  # ub_scb <- psi + qtmax * sqrt(psi_var)

  return(list(Estimates = psi,
              Variances = psi_var,
              CI_LB = lb,
              CI_UB = ub)) #,
              # SCB_LB = lb_scb,
              # SCB_UB = ub_scb))
}


