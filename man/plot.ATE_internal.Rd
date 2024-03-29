% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3methods.R
\name{plot.ATE_internal}
\alias{plot.ATE_internal}
\title{Plot method for objects of class "ATE_internal"}
\usage{
\method{plot}{ATE_internal}(x, source_names, ...)
}
\arguments{
\item{x}{Object of class "ATE_internal".}

\item{source_names}{optional, vector of character strings specifying the names of the sources. Defaults are the values in \code{S} provided by the user to \code{\link{ATE_internal}}.}

\item{...}{Other arguments, which are passed to \code{\link[metafor]{forest.rma}}.}
}
\value{
No value is returned.
}
\description{
This function creates forest plots of objects of class "ATE_internal".
}
\examples{
\donttest{
ai <- ATE_internal(
  X = dat_multisource[, 1:10],
  Y = dat_multisource$Y,
  S = dat_multisource$S,
  A = dat_multisource$A,
  source_model = "MN.glmnet",
  source_model_args = list(),
  treatment_model_type = "separate",
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
plot(ai)
}

}
\seealso{
\code{\link{ATE_internal}}
}
