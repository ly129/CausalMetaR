#' Plot method for objects of class "STE_int"
#'
#' This function creates forest plots of objects of class "STE_int".
#'
#' @param x Object of class "STE_int".
#' @param header character string specifying the title of the column of subgroup names
#' @param xlab character string specifying the title for the x-axis.
#' @param include_scb logical scalar specifying whether to include simultaneous 95\% confidence bands in the forest plot
#' @param ... Other arguments, which are passed to \code{\link[metafor]{forest.rma}}.
#' @return No value is returned.
#' @seealso \code{\link{STE_int}}
#'
#'
#'
#' @export

plot.STE_int <- function(x, header = 'Subgroup', xlab = 'Treatment Effect',
                       include_scb = TRUE, ...){
  if (!inherits(x, "STE_int")){
    stop("Argument 'x' must be an object of class \"STE_int\".")
  }

  forest(x = x$plot_psi,
         vi = x$plot_psi_var,
         slab = x$snames,
         ilab = x$xtildenames,
         ilab.xpos = -3,
         header = header,
         xlab = xlab,
         ...)

  if (include_scb){
    for (i in 1:(x$n_x_tilde * x$no_S)) {
      graphics::segments(x0 = x$plot_scb[i, 1],
                         y0 = x$n_x_tilde * x$no_S + 1 - i,
                         x1 = x$plot_scb[i, 2])
    }
  }
}




#' Print method for objects of class "STE_int"
#'
#' Print method for objects of class "STE_int"
#'
#' @param x Object of class "STE_int".
#' @param digits Integer specifying the number of decimal places to display.
#' @param ... Other arguments.
#' @return No value is returned.
#' @seealso \code{\link{STE_int}}
#'
#'
#'
#' @export

print.STE_int <- function(x, digits = 4, ...){
  if (!inherits(x, "STE_int")){
    stop("Argument 'x' must be an object of class \"STE_int\".")
  }
  cat('SUBGROUP TREATMENT EFFECT ESTIMATES')
  cat("\n-----------------------------------\n\n")

  df <- data.frame(Study = rep(1:x$no_S, each = x$n_x_tilde),
                   Subgroup = rep(1:x$n_x_tilde, times = x$no_S),
                   Estimate = x$plot_psi,
                   SE = sqrt(x$plot_psi_var),
                   ci.lb = x$plot_psi - qnorm(p = 0.975) * sqrt(x$plot_psi_var),
                   ci.ub = x$plot_psi + qnorm(p = 0.975) * sqrt(x$plot_psi_var))
  colnames(df)[4:6] <- c('Std. Error', 'Lower CI', 'Upper CI')
  print(df, row.names = FALSE, digits = digits)
}
