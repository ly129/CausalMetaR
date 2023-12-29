#' Multi-source dataset
#'
#' Simulated multi-source dataset consisting of 3 sources.
#'
#' @docType data
#'
#' @format A data frame with 3,917 rows and 13 columns. The columns are:
#' \tabular{ll}{
#'  \code{EM} \tab Effect modifier. \cr
#'  \code{X2}, ..., \code{X10} \tab Covariates. \cr
#'  \code{S} \tab Source indicator. \cr
#'  \code{A} \tab Treatment (1 for treated and 0 for untreated). \cr
#'  \code{Y} \tab Outcome.
#' }
"dat_multisource"



#' External dataset
#'
#' Simulated external dataset.
#'
#' @docType data
#'
#' @format A data frame with 10,083 rows and 13 columns. The columns are:
#' \tabular{ll}{
#'  \code{EM} \tab Effect modifier. \cr
#'  \code{X2}, ..., \code{X10}  \tab  Covariates. \cr
#'  \code{S} \tab Source indicator. \cr
#'  \code{A} \tab Treatment (1 for treated and 0 for untreated). \cr
#'  \code{Y} \tab Outcome.
#' }
"dat_external"
