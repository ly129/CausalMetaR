#' Nested data
#'
#' Simulated nested dataset for CMetafoR
#'
#' @format `dat_nested`
#' A data frame with 1019 rows and 16 columns
#' \describe{
#'  \item{EM}{Effect modifier}
#'  \item{X2, ..., X10}{9 covariates}
#'  \item{R}{Indicator for nested data: R = 1 for nested data (not used)}
#'  \item{S}{Source}
#'  \item{A}{Treatment: 1 - treated; 0 - untreated}
#'  \item{Y}{Outcome}
#'  \item{Y1, Y0}{Counterfactual outcomes (not used)}
#' }
"dat_nested"



#' External data
#'
#' Simulated external dataset for CMetafoR
#'
#' @format `dat_external`
#' A data frame with 8981 rows and 16 columns
#' \describe{
#'  \item{EM}{Effect modifier}
#'  \item{X2, ..., X10}{9 covariates}
#'  \item{R}{Indicator for nested data: R = 0 for external data (not used)}
#'  \item{S}{Source}
#'  \item{A}{Treatment: 1 - treated; 0 - untreated}
#'  \item{Y}{Outcome}
#'  \item{Y1, Y0}{Counterfactual outcomes (not used)}
#' }
"dat_external"
