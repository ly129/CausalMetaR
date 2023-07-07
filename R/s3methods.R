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

  forest(x = x$df_dif$Estimate,
         vi = x$df_dif$SE^2,
         slab = x$snames,
         ilab = x$xtildenames,
         ilab.xpos = -3,
         header = header,
         xlab = xlab,
         ...)

  if (include_scb){
    no_S <- length(unique(x$df_dif$Study))
    no_x_tilde <- length(unique(x$df_dif$Subgroup))

    for (i in 1:(no_x_tilde * no_S)) {
      graphics::segments(x0 = x$df_dif$scb.lb[i],
                         y0 = no_x_tilde * no_S + 1 - i,
                         x1 = x$df_dif$scb.ub[i])
    }
  }
}




#' Print method for objects of class "ATE_int", "ATE_ext", "STE_int", or "STE_ext"
#'
#' Print method for objects of class "ATE_int", "ATE_ext", "STE_int", or "STE_ext"
#'
#' @param x Object of class "ATE_int", "ATE_ext", "STE_int", or "STE_ext".
#' @param digits Integer specifying the number of decimal places to display.
#' @param ... Other arguments (ignored).
#' @return No value is returned.
#' @seealso \code{\link{ATE_int}}, \code{\link{ATE_ext}}, \code{\link{STE_int}}, \code{\link{STE_ext}}
#'
#'
#'
#' @export

print.STE_int <- function(x, digits = 4, ...){
  if (!inherits(x, "STE_int")){
    stop("Argument 'x' must be an object of class \"STE_int\".")
  }

  no_S <- length(unique(x$df_dif$Study))
  no_x_tilde <- length(unique(x$df_dif$Subgroup))

  study_num <- matrix(1:no_S, nrow = 1)
  append_row <- matrix(NA, ncol = no_S, nrow = no_x_tilde - 1)
  study_lab <- c(rbind(study_num, append_row))

  df_dif <- x$df_dif
  df_dif$Study <- study_lab

  cat('SUBGROUP TREATMENT EFFECT ESTIMATES IN INTERNAL POPULATIONS\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(df_dif, digits = digits, ATE = FALSE, internal = TRUE)
}

#' @rdname print.STE_int
#' @export
print.ATE_int <- function(x, digits = 4, ...){
  if (!inherits(x, "ATE_int")){
    stop("Argument 'x' must be an object of class \"ATE_int\".")
  }

  cat('AVERAGE TREATMENT EFFECT ESTIMATES IN INTERNAL POPULATIONS\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(x$df_dif, digits = digits, ATE = TRUE, internal = TRUE)
}


#' @rdname print.STE_int
#' @export
print.STE_ext <- function(x, digits = 4, ...){
  if (!inherits(x, "STE_ext")){
    stop("Argument 'x' must be an object of class \"STE_ext\".")
  }

  cat('SUBGROUP TREATMENT EFFECT ESTIMATES IN AN EXTERNAL POPULATION\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(x$df_dif, digits = digits, ATE = FALSE, internal = FALSE)
}

#' @rdname print.STE_int
#' @export
print.ATE_ext <- function(x, digits = 4, ...){
  if (!inherits(x, "ATE_ext")){
    stop("Argument 'x' must be an object of class \"ATE_ext\".")
  }

  cat('AVERAGE TREATMENT EFFECT ESTIMATES IN AN EXTERNAL POPULATION\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(x$df_dif, digits = digits, ATE = TRUE, internal = FALSE)
}

#' Summary method for objects of class "ATE_int", "ATE_ext", "STE_int", or "STE_ext"
#'
#' Summary method for objects of class "ATE_int", "ATE_ext", "STE_int", or "STE_ext"
#'
#' @param object Object of class "ATE_int", "ATE_ext", "STE_int", or "STE_ext".
#' @param digits Integer specifying the number of decimal places to display.
#' @param ... Other arguments.
#' @return No value is returned.
#' @seealso \code{\link{ATE_int}}, \code{\link{ATE_ext}}, \code{\link{STE_int}}, \code{\link{STE_ext}}
#'
#'
#'
#' @export

summary.STE_int <- function(object, digits = 4, ...){
  if (!inherits(object, "STE_int")){
    stop("Argument 'object' must be an object of class \"STE_int\".")
  }

  no_S <- length(unique(object$df_dif$Study))
  no_x_tilde <- length(unique(object$df_dif$Subgroup))

  study_num <- matrix(1:no_S, nrow = 1)
  append_row <- matrix(NA, ncol = no_S, nrow = no_x_tilde - 1)
  study_lab <- c(rbind(study_num, append_row))

  df_dif <- object$df_dif; df_dif$Study <- study_lab
  df_A0 <- object$df_A0; df_A0$Study <- study_lab
  df_A1 <- object$df_A1; df_A1$Study <- study_lab

  cat('SUBGROUP TREATMENT EFFECT ESTIMATES IN INTERNAL POPULATIONS\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(df_dif, digits = digits, ATE = FALSE, internal = TRUE)

  cat('\n\nPotential outcome mean estimates under A = 0:\n')
  cat("---------------------------------------------\n")
  my_print(df_A0, digits = digits, ATE = FALSE, internal = TRUE)

  cat('\n\nPotential outcome mean estimates under A = 1:\n')
  cat("---------------------------------------------\n")
  my_print(df_A1, digits = digits, ATE = FALSE, internal = TRUE)
}


#' @rdname summary.STE_int
#' @export
summary.STE_ext <- function(object, digits = 4, ...){
  if (!inherits(object, "STE_ext")){
    stop("Argument 'object' must be an object of class \"STE_ext\".")
  }

  cat('SUBGROUP TREATMENT EFFECT ESTIMATES IN AN EXTERNAL POPULATION\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(object$df_dif, digits = digits, ATE = FALSE, internal = FALSE)

  cat('\n\nPotential outcome mean estimates under A = 0:\n')
  cat("---------------------------------------------\n")
  my_print(object$df_A0, digits = digits, ATE = FALSE, internal = FALSE)

  cat('\n\nPotential outcome mean estimates under A = 1:\n')
  cat("---------------------------------------------\n")
  my_print(object$df_A1, digits = digits, ATE = FALSE, internal = FALSE)
}

#' @rdname summary.STE_int
#' @export
summary.ATE_ext <- function(object, digits = 4, ...){
  if (!inherits(object, "ATE_ext")){
    stop("Argument 'object' must be an object of class \"ATE_ext\".")
  }

  cat('AVERAGE TREATMENT EFFECT ESTIMATES IN AN EXTERNAL POPULATION\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(object$df_dif, digits = digits, ATE = TRUE, internal = FALSE)

  cat('\n\nPotential outcome mean estimates under A = 0:\n')
  cat("---------------------------------------------\n")
  my_print(object$df_A0, digits = digits, ATE = TRUE, internal = FALSE)

  cat('\n\nPotential outcome mean estimates under A = 1:\n')
  cat("---------------------------------------------\n")
  my_print(object$df_A1, digits = digits, ATE = TRUE, internal = FALSE)
}

#' @rdname summary.STE_int
#' @export
summary.ATE_int <- function(object, digits = 4, ...){
  if (!inherits(object, "ATE_int")){
    stop("Argument 'object' must be an object of class \"ATE_int\".")
  }

  cat('AVERAGE TREATMENT EFFECT ESTIMATES IN AN INTERNAL POPULATIONS\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(object$df_dif, digits = digits, ATE = TRUE, internal = TRUE)

  cat('\n\nPotential outcome mean estimates under A = 0:\n')
  cat("---------------------------------------------\n")
  my_print(object$df_A0, digits = digits, ATE = TRUE, internal = TRUE)

  cat('\n\nPotential outcome mean estimates under A = 1:\n')
  cat("---------------------------------------------\n")
  my_print(object$df_A1, digits = digits, ATE = TRUE, internal = TRUE)
}



my_print <- function(df, digits, ATE, internal, ...){
  my_fun <- function(x){
    sprintf(paste0("%.", digits, "f"), round(x, digits))
  }
  mycols <- which(colnames(df) %in% c('Study', 'Subgroup'))

  if (ATE & !internal){
    df_num <- df
  } else if (ATE & internal){
    df_num <- df[, -mycols]
    df_char <- data.frame(Study = df[, mycols])
    df_char <- format.data.frame(data.frame(lapply(df_char, as.character)), na.encode = FALSE)
  } else if (!ATE & !internal){
    df_num <- df[, -mycols]
    df_char <- data.frame(Subgroup = df[, mycols])
    df_char <- format.data.frame(data.frame(lapply(df_char, as.character)), na.encode = FALSE)
  } else if (!ATE & internal){
    df_num <- df[, -mycols]
    df_char <- df[, mycols]
    df_char <- format.data.frame(data.frame(lapply(df_char, as.character)), na.encode = FALSE)
  }

  df_rounded <- format.data.frame(data.frame(lapply(df_num, my_fun)))
  if (!ATE){
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
