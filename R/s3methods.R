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
    for (i in 1:(x$no_x_tilde * x$no_S)) {
      graphics::segments(x0 = x$plot_scb[i, 1],
                         y0 = x$no_x_tilde * x$no_S + 1 - i,
                         x1 = x$plot_scb[i, 2])
    }
  }
}




#' Print method for objects of class "STE_int" or "STE_ext"
#'
#' Print method for objects of class "STE_int" or "STE_ext"
#'
#' @param x Object of class "STE_int" or "STE_ext".
#' @param digits Integer specifying the number of decimal places to display.
#' @param ... Other arguments (ignored).
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

  study_num <- matrix(1:x$no_S, nrow = 1)
  append_row <- matrix(NA, ncol = x$no_S, nrow = x$no_x_tilde - 1)
  study_lab <- c(rbind(study_num, append_row))

  df <- data.frame(Study = study_lab,
                   Subgroup = rep(1:x$no_x_tilde, times = x$no_S),
                   Estimate = x$plot_psi,
                   SE = sqrt(x$plot_psi_var),
                   ci.lb = x$plot_psi - qnorm(p = 0.975) * sqrt(x$plot_psi_var),
                   ci.ub = x$plot_psi + qnorm(p = 0.975) * sqrt(x$plot_psi_var),
                   scb.lb = x$plot_scb[, 1],
                   scb.ub = x$plot_scb[, 2])

  cat('SUBGROUP TREATMENT EFFECT ESTIMATES IN INTERNAL POPULATIONS\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(df, digits = digits, include_scb = TRUE)
}

#' @rdname print.STE_int
#' @export
print.ATE_int <- function(x, digits = 4, ...){
  if (!inherits(x, "ATE_int")){
    stop("Argument 'x' must be an object of class \"ATE_int\".")
  }

  df <- data.frame(Study = 1:x$no_S,
                   Estimate = unname(x$plot_psi),
                   SE = unname(sqrt(x$plot_psi_var)),
                   ci.lb = unname(x$plot_psi_CI[, 1]),
                   ci.ub = unname(x$plot_psi_CI[, 2]))

  cat('AVERAGE TREATMENT EFFECT ESTIMATES IN INTERNAL POPULATIONS\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(df, digits = digits, include_scb = FALSE)
}


#' @rdname print.STE_int
#' @export
print.STE_ext <- function(x, digits = 4, ...){
  if (!inherits(x, "STE_ext")){
    stop("Argument 'x' must be an object of class \"STE_ext\".")
  }

  df <- data.frame(Estimate = x$plot_phi,
                   SE = sqrt(x$plot_phi_var),
                   ci.lb = x$plot_phi_CI[1],
                   ci.ub = x$plot_phi_CI[2])

  cat('SUBGROUP TREATMENT EFFECT ESTIMATES IN AN EXTERNAL POPULATION\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(df, digits = digits, include_scb = FALSE)
}

#' @rdname print.STE_int
#' @export
print.ATE_ext <- function(x, digits = 4, ...){
  if (!inherits(x, "ATE_ext")){
    stop("Argument 'x' must be an object of class \"ATE_ext\".")
  }

  df <- data.frame(Estimate = x$plot_phi,
                   SE = sqrt(x$plot_phi_var),
                   ci.lb = x$plot_phi_CI[1],
                   ci.ub = x$plot_phi_CI[2])

  cat('AVERAGE TREATMENT EFFECT ESTIMATES IN AN EXTERNAL POPULATION\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(df, digits = digits, include_scb = FALSE)
}

#' Summary method for objects of class "STE_int"
#'
#' Summary method for objects of class "STE_int"
#'
#' @param object Object of class "STE_int".
#' @param digits Integer specifying the number of decimal places to display.
#' @param ... Other arguments.
#' @return No value is returned.
#' @seealso \code{\link{STE_int}}
#'
#'
#'
#' @export

summary.STE_int <- function(object, digits = 4, ...){
  if (!inherits(object, "STE_int")){
    stop("Argument 'object' must be an object of class \"STE_int\".")
  }

  study_num <- matrix(1:object$no_S, nrow = 1)
  append_row <- matrix(NA, ncol = object$no_S, nrow = object$no_x_tilde - 1)
  study_lab <- c(rbind(study_num, append_row))

  df <- data.frame(Study = study_lab,
                   Subgroup = rep(1:object$no_x_tilde, times = object$no_S),
                   Estimate = object$plot_psi,
                   SE = sqrt(object$plot_psi_var),
                   ci.lb = object$plot_psi - qnorm(p = 0.975) * sqrt(object$plot_psi_var),
                   ci.ub = object$plot_psi + qnorm(p = 0.975) * sqrt(object$plot_psi_var),
                   scb.lb = object$plot_scb[, 1],
                   scb.ub = object$plot_scb[, 2])

  df_A1 <- df_A0 <- data.frame(Study = study_lab,
                               Subgroup = rep(1:object$no_x_tilde, times = object$no_S),
                               Estimate = NA,
                               SE = NA,
                               ci.lb = NA,
                               ci.ub = NA,
                               scb.lb = NA,
                               scb.ub = NA)
  row_ind <- 1
  for (i in 1:object$no_S){
    for (j in 1:object$no_x_tilde){
      df_A1[row_ind, 'Estimate'] <- object[[i]]$Estimates[j, 1]
      df_A0[row_ind, 'Estimate'] <- object[[i]]$Estimates[j, 2]

      df_A1[row_ind, 'SE'] <- sqrt(object[[i]]$Variances[j, 1])
      df_A0[row_ind, 'SE'] <- sqrt(object[[i]]$Variances[j, 2])

      df_A1[row_ind, 'ci.lb'] <- object[[i]]$CI_LB[j, 1]
      df_A0[row_ind, 'ci.lb'] <- object[[i]]$CI_LB[j, 2]

      df_A1[row_ind, 'ci.ub'] <- object[[i]]$CI_UB[j, 1]
      df_A0[row_ind, 'ci.ub'] <- object[[i]]$CI_UB[j, 2]

      df_A1[row_ind, 'scb.lb'] <- object[[i]]$SCB_LB[j, 1]
      df_A0[row_ind, 'scb.lb'] <- object[[i]]$SCB_LB[j, 2]

      df_A1[row_ind, 'scb.ub'] <- object[[i]]$SCB_UB[j, 1]
      df_A0[row_ind, 'scb.ub'] <- object[[i]]$SCB_UB[j, 2]

      row_ind <- row_ind + 1
    }
  }

  cat('SUBGROUP TREATMENT EFFECT ESTIMATES IN INTERNAL POPULATIONS\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(df, digits = digits, include_scb = TRUE)

  cat('\n\nPotential outcome mean estimates under A = 0:\n')
  cat("-------------------------------------------\n")
  my_print(df_A0, digits = digits, include_scb = TRUE)

  cat('\n\nPotential outcome mean estimates under A = 1:\n')
  cat("-------------------------------------------\n")
  my_print(df_A1, digits = digits, include_scb = TRUE)
}



my_print <- function(df, digits, include_scb, ...){
  my_fun <- function(x){
    sprintf(paste0("%.", digits, "f"), x)
  }
  mycols <- which(colnames(df) %in% c('Study', 'Subgroup'))

  if (length(mycols) > 0){
    df_rounded <- round(df[, -mycols], digits = digits)
    df_char <- df[, mycols]
    df_char <- format.data.frame(data.frame(lapply(df_char, as.character)), na.encode = FALSE)
  } else {
    df_rounded <- round(df, digits = digits)
  }

  df_rounded <- format.data.frame(data.frame(lapply(df_rounded, my_fun)))
  if (include_scb){
    colnames(df_rounded)[2:6] <- c('SE', 'Lower 95% CI', 'Upper 95% CI', 'Lower 95% SCB', 'Upper 95% SCB')
  } else {
    colnames(df_rounded)[2:4] <- c('SE', 'Lower 95% CI', 'Upper 95% CI')
  }

  if (length(mycols) > 0){
    out <- cbind(df_char, df_rounded)
  } else {
    out <- df_rounded
  }
  print(out, na.print = "", row.names = FALSE)
}
