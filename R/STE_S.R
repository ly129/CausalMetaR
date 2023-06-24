#' Transporting subgroup treatment effects from multi-source population to an internal source-specific population
#'
#' @description
#' Doubly-robust and efficient estimator for the average treatment effect of an internal source-specific target population using \eqn{m} multi-source data.   
#' 
#' @param X The covariate matrix/data frame with \eqn{n=n_1+...+n_m} rows and q coloums.
#' @param Y The (binary/categorical/continuous) outcome, which is a length \eqn{n} vector.
#' @param S The (numeric) source which is a length \eqn{n} vector. 
#' @param A The (binary) treatment, which is a length \eqn{n} vector.
#' @param source_model The multi-nomial model for estimating \eqn{P(S=s|X)}. It has two options: \code{SL.glmnet.multinom} and \code{SL.nnet.multinom}. The default is \code{SL.glmnet.multinom}.
#' @param source_model_args The arguments (in \pkg{SuperLearner}) for the source model.
#' @param treatment_model_type The options for how the treatment_model \eqn{P(A=1|X, S=s)} is estimated. It includes \code{separate} and \code{joint}, with the default being \code{separate}. When \code{separate} is selected, 
#' \eqn{P(A=1|X, S=s)} is estimated by fitting the model (regressing \eqn{A} on \eqn{X}) within each specific internal source population (S=s). When \code{joint} is selected, \eqn{P(A=1|X, S=s)} 
#' is estimated by fitting the model (regressing \eqn{A} on \eqn{X} and \eqn{S}) using the multi-source population and then estimating the probability by fitting the model while suppressing the \eqn{S=s}. 
#' In both cases, the propensity score is calculated as $P(A=1|X)=\sum_{s=1}^{m}P(A=1|X, S=s)P(S=s|X)$.
#' @param treatment_model The treatment model \eqn{P(A=1|X, S=s)} is estimated using \pkg{SuperLearner}. If, for example, the preference is to use only logistic regression for the probability estimation, 
#' please ensure that only \code{glm} is included in the \pkg{SuperLearner} library within the \code{treatment_model_args}.
#' @param treatment_model_args The arguments (in \pkg{SuperLearner}) for the treatment model.
#' @param outcome_model The same as \code{treatment_model}.
#' @param outcome_model_args The arguments (in \pkg{SuperLearner}) for the outcome model.
#'
#' @details
#' Data structure: multi-source data contain outcome Y, source S, treatment A, and covariates X ($n \times p$); external data contain only covariate X0 ($n0 \times p$).
#' The estimator is doubly robust and non-parametrically efficient. Three nuisance parameters are estimated, 
#' the source model $q_s(X)=P(S=s|X)$, the propensity score model $\eta_a(X)=P(A=a|X)$, and the outcome model $\mu_a(X)=E(Y|X, A=a)$. The nuisance parameters are allowed to be estimated by \pkg{SuperLearner}. The estimator is
#' $$
#'  \dfrac{\widehat \kappa}{n}\sum\limits_{i=1}^{n} \Bigg[ I(S_i = s) \widehat \mu_a(X_i) 
#'  +I(A_i = a) \dfrac{\widehat q_{s}(X_i)}{\widehat \eta_a(X_i)}  \Big\{ Y_i - \widehat \mu_a(X_i) \Big\} \Bigg],
#' $$
#' where $\widehat \kappa=\{n^{-1} \sum_{i=1}^n I(S_i=s)\}^{-1}$.
#' To achieve the non-parametrical efficiency and asymptotic normality, it requires that $||\widehat \mu_a(X) -\mu_a(X)||\big\{||\widehat \eta_a(X) -\eta_a(X)||+||\widehat q_s(X) -q_s(X)||\big\}=o_p(n^{-1/2})$. 
#' In addition, to avoid the Donsker class assumption, the estimation is done by sample splitting and cross-fitting.
#' When one source of data is a randomized trial, it is still recommended to estimate the propensity score for optimal efficiency. 
#' Since the non-parametric influence function is the same as the efficient semi-parametric efficient influence function when the propensity score is known and incorporating the assumption $Y\prep S|(X, A=a)$, the inference stays the same. 
#'
#' @return A list with the following four elements.
#'   \item{Estimates}{The point estimate of the ATE for each of s.}
#'   \item{Variances}{The asymptotic variances of the point estimates, which are calculated based on the (efficient) influence function.}
#'   \item{CI_LB}{The lower bounds of the 95% confidence intervals.}
#'   \item{CI_UB}{The upper bounds of the 95% confidence intervals.}
#'
#' @references Dahabreh, I.J., Robertson, S.E., Petito, L.C., Hernán, M.A. and Steingrimsson, J.A.. (2019) \emph{Efficient and robust methods for causally 
#' interpretable meta‐analysis: Transporting inferences from multiple randomized trials to a target population}, Biometrics.
#' 
#' @examples 
#'
CMetafoR.STE.S <- function(
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

  unique_X <- sort(unique(X[, 1]))
  n_x_tilde <- length(unique_X)
  output <- vector(mode = "list", length = n_x_tilde)
  names(output) <- paste(names(X)[1], "=", unique_X)

  plot_psi <- plot_psi_var <- numeric(n_x_tilde * no_S)
  plot_scb <- matrix(nrow = n_x_tilde * no_S, ncol = 2)

  for (i in 1:n_x_tilde) {
    x_tilde <- unique_X[i]
    for (s in unique_S) {
      tmp1 <- tmp2 <- matrix(0, nrow = n, ncol = 2)
      I_xs <- which((X[, 1] == x_tilde) & (S == s))
      kappa <- 1/(length(I_xs)/n)
      tmp1[I_xs, 1] <- pred_Y1[I_xs]
      tmp1[I_xs, 2] <- pred_Y0[I_xs]

      I_xa1 <- which((X[, 1] == x_tilde) & (A == 1))
      I_xa0 <- which((X[, 1] == x_tilde) & (A == 0))

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

    tmax <- apply(abs(matrix(rnorm(length(unique(X[,1])) * 1e6),
                             nrow = length(unique(X[,1])), ncol = 1e6)), 2, max)
    qtmax <- quantile(tmax, 0.95)

    lb_scb <- psi - qtmax * sqrt(psi_var)
    ub_scb <- psi + qtmax * sqrt(psi_var)

    output[[i]] <- list(Estimates = psi,
                        Variances = psi_var,
                        CI_LB = lb,
                        CI_UB = ub,
                        SCB_LB = lb_scb,
                        SCB_UB = ub_scb)
    plot_psi[((i - 1) * no_S + 1):(i * no_S)] <- psi[, 1] - psi[, 2]
    plot_psi_var[((i - 1) * no_S + 1):(i * no_S)] <- psi_var[,1] + psi_var[,2]
    plot_scb[((i - 1) * no_S + 1):(i * no_S), 1] <- psi[, 1] - psi[, 2] - qtmax * sqrt(psi_var[,1] + psi_var[,2])
    plot_scb[((i - 1) * no_S + 1):(i * no_S), 2] <- psi[, 1] - psi[, 2] + qtmax * sqrt(psi_var[,1] + psi_var[,2])
  }

  snames <- rep(paste("Study =", unique_S), n_x_tilde)
  xtildenames <- character(length = n_x_tilde * no_S)
  xtildenames[1:(n_x_tilde * no_S) %% no_S == 1] <- c(paste(names(X)[1], "=", unique_X))

  forest(x = plot_psi,
         vi = plot_psi_var,
         slab = xtildenames,
         ilab = snames,
         ilab.xpos = -3,
         header = "Subgroup",
         xlab = "Treatment Effect")
  for (i in 1:(n_x_tilde * no_S)) {
    segments(x0 = plot_scb[i, 1],
             y0 = n_x_tilde * no_S + 1 - i,
             x1 = plot_scb[i, 2])
  }

  return(output)
}


