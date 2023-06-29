#' Plot method for objects of class "STE_s"
#'
#' This function creates forest plots of objects of class "STE_S".
#'
#' @param x Object of class "STE_S".
#' @param header character string specifying the title of the column of subgroup names
#' @param xlab character string specifying the title for the x-axis.
#' @param ... Other arguments, which are passed to \code{\link[metafor]{forest.rma}}.
#' @return No value is returned.
#' @seealso \code{\link{CMetafoR.STE.S}}
#'
#'
#'
#' @export

plot.STE_S <- function(x, header = 'Subgroup', xlab = 'Treatment Effect', ...){
  if (!inherits(x, "STE_S")){
    stop("Argument 'x' must be an object of class \"STE_S\".")
  }

  forest(x = x$plot_psi,
         vi = x$plot_psi_var,
         slab = x$snames,
         ilab = x$xtildenames,
         ilab.xpos = -3,
         header = header,
         xlab = xlab,
         ...)

  for (i in 1:(x$n_x_tilde * x$no_S)) {
    graphics::segments(x0 = x$plot_scb[i, 1],
                       y0 = x$n_x_tilde * x$no_S + 1 - i,
                       x1 = x$plot_scb[i, 2])
  }
}
