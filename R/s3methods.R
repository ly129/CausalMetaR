#' Plot method for objects of class "STE_internal"
#'
#' This function creates forest plots of objects of class "STE_internal".
#'
#' Note that users may need to custom set the argument \code{ilab.xpos} which specifies the position (along the x-axis) of the effect modifier header and subgroup labels. See \code{\link[metafor]{forest.rma}} for further details.
#'
#' @param x Object of class "STE_internal".
#' @param use_scb logical scalar specifying whether the intervals in the forest plot should be simultaneous confidence bands (rather than confidence intervals). The default is \code{FALSE}.
#' @param header optional, vector of character strings of length 3, headers for the source, effect modifier subgroup and the estimates in the forest plot.
#' @param source_names optional, vector of character strings specifying the names of the sources. Defaults are the values in \code{S} provided by the user to \code{\link{STE_internal}}.
#' @param subgroup_names optional, vector of character strings specifying the names of the effect modifier subgroups. Defaults are the values in \code{EM} provided by the user to \code{\link{STE_internal}}.
#' @param ... Other arguments, which are passed to \code{\link[metafor]{forest.rma}}.
#' @return No value is returned.
#' @seealso \code{\link{STE_internal}}
#'
#' @examples
#' \donttest{
#' si <- STE_internal(
#'   X = dat_multisource[, 2:10],
#'   Y = dat_multisource$Y,
#'   EM = dat_multisource$EM,
#'   S = dat_multisource$S,
#'   A = dat_multisource$A,
#'   cross_fitting = FALSE,
#'   source_model = "MN.nnet",
#'   source_model_args = list(),
#'   treatment_model_type = "separate",
#'   treatment_model_args = list(
#'     family = binomial(),
#'     SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
#'     cvControl = list(V = 5L)
#'   ),
#'   outcome_model_args = list(
#'     family = gaussian(),
#'     SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
#'     cvControl = list(V = 5L)
#'   )
#' )
#' plot(si)
#' }
#'
#' @export

plot.STE_internal <- function(x,
                            use_scb = FALSE,
                            header = c("Source",
                                       "Subgroup",
                                       ifelse(use_scb,
                                              "Estimate [95% SCB]",
                                              "Estimate [95% CI]")),
                            source_names,
                            subgroup_names,
                            ...){
  if (!inherits(x, "STE_internal")){
    stop("Argument 'x' must be an object of class \"STE_internal\".")
  }
  all_args <- as.list(match.call())[-1]
  args <- all_args[!names(all_args) %in% c('x', 'use_scb', 'source_names', 'subgroup_names')]

  no_S <- length(x$source_names)
  no_EM <- length(x$subgroup_names)

  slab <- character(length = no_EM * no_S)
  if (missing(source_names)){
    slab[1:(no_EM * no_S) %% no_EM == 1] <- x$source_names
    args$slab <- slab
  } else {
    if (length(source_names) != no_S){
      stop(paste0('The length of source_names does not match the number of sources (', no_S, ')'))
    }
    slab[1:(no_EM * no_S) %% no_EM == 1] <- source_names
    args$slab <- slab
  }

  if (missing(subgroup_names)){
    args$ilab <- rep(x$subgroup_names, times = no_S)
  } else {
    if (length(subgroup_names) != no_EM){
      stop(paste0('The length of subgroup_names does not match the number of subgroups (', no_EM, ')'))
    }
    args$ilab <- rep(subgroup_names, times = no_S)
  }
  if (!('shade' %in% names(args))){
    shade_temp <- c(rep(TRUE, times = no_EM), rep(FALSE, times = no_EM))
    args$shade <- rep(shade_temp, length.out = no_S * no_EM)
  }
  if (!('xlab' %in% names(args))) {
    args$xlab <- 'Treatment Effect'
  }
  if (!('refline' %in% names(args))){
    args$refline <- NA
  }
  args$x <- x$df_dif$Estimate

  if (!use_scb){
    args$vi <- x$df_dif$SE^2
  } else {
    args$vi <- ((x$df_dif$scb.ub - x$df_dif$scb.lb) / (2 * qnorm(0.975)))^2
    if ('level' %in% names(args)){
      if (args$level != 0.95){
        stop("The argument 'level' cannot be set to a value other than 0.95 when using simultaneous confidence bands.")
      }
    }
  }

  args$header <- c(header[1], header[3])

  if (!('ilab.xpos' %in% names(args))){
    min_lb <- min(args$x - qnorm(0.975) * sqrt(args$vi))
    median_est <- median(args$x)
    args$ilab.xpos <- min_lb - (median_est - min_lb) # / 2
  }

  do.call(forest, args)

  text(x = args$ilab.xpos,
       y = no_EM * no_S + 2.01,
       labels = header[2],
       font = 2)
}


#' Plot method for objects of class "ATE_internal"
#'
#' This function creates forest plots of objects of class "ATE_internal".
#'
#' @param x Object of class "ATE_internal".
#' @param source_names optional, vector of character strings specifying the names of the sources. Defaults are the values in \code{S} provided by the user to \code{\link{ATE_internal}}.
#' @param ... Other arguments, which are passed to \code{\link[metafor]{forest.rma}}.
#' @return No value is returned.
#' @seealso \code{\link{ATE_internal}}
#'
#' @examples
#' \donttest{
#' ai <- ATE_internal(
#'   X = dat_multisource[, 1:10],
#'   Y = dat_multisource$Y,
#'   S = dat_multisource$S,
#'   A = dat_multisource$A,
#'   source_model = "MN.glmnet",
#'   source_model_args = list(),
#'   treatment_model_type = "separate",
#'   treatment_model_args = list(
#'     family = binomial(),
#'     SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
#'     cvControl = list(V = 5L)
#'   ),
#'   outcome_model_args = list(
#'     family = gaussian(),
#'     SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
#'     cvControl = list(V = 5L)
#'   )
#' )
#' plot(ai)
#' }
#'
#' @export

plot.ATE_internal <- function(x,
                            source_names,
                            ...){
  if (!inherits(x, "ATE_internal")){
    stop("Argument 'x' must be an object of class \"ATE_internal\".")
  }

  all_args = as.list(match.call())[-1]
  args <- all_args[!names(all_args) %in% c('x')]

  no_S <- length(unique(x$df_dif$Source))

  if (missing(source_names)){
    args$slab <- x$source_names
  } else {
    if (length(source_names) != no_S){
      stop(paste0('The length of source_names does not match the number of sources (', no_S, ')'))
    }
    args$slab <- source_names
  }

  if (!('shade' %in% names(args))){
    args$shade <- 'zebra'
  }
  if (!('xlab' %in% names(args))) {
    args$xlab <- 'Treatment Effect'
  }
  if (!('refline' %in% names(args))){
    args$refline <- NA
  }
  if (!('header' %in% names(args))){
    args$header <- 'Source'
  }
  args$x <- x$df_dif$Estimate
  args$vi <- x$df_dif$SE^2
  do.call(forest, args)
}




#' Print method for objects of class "ATE_internal", "ATE_external", "STE_internal", or "STE_external"
#'
#' Print method for objects of class "ATE_internal", "ATE_external", "STE_internal", or "STE_external"
#'
#' @param x Object of class "ATE_internal", "ATE_external", "STE_internal", or "STE_external".
#' @param digits Integer specifying the number of decimal places to display.
#' @param ... Other arguments (ignored).
#' @return No value is returned.
#' @seealso \code{\link{ATE_internal}}, \code{\link{ATE_external}}, \code{\link{STE_internal}}, \code{\link{STE_external}}
#'
#' @examples
#' \donttest{
#' si <- STE_internal(
#'   X = dat_multisource[, 2:10],
#'   Y = dat_multisource$Y,
#'   EM = dat_multisource$EM,
#'   S = dat_multisource$S,
#'   A = dat_multisource$A,
#'   cross_fitting = FALSE,
#'   source_model = "MN.nnet",
#'   source_model_args = list(),
#'   treatment_model_type = "separate",
#'   treatment_model_args = list(
#'     family = binomial(),
#'     SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
#'     cvControl = list(V = 5L)
#'   ),
#'   outcome_model_args = list(
#'     family = gaussian(),
#'     SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
#'     cvControl = list(V = 5L)
#'   )
#' )
#' print(si)
#' }
#'
#' @export

print.STE_internal <- function(x, digits = 4, ...){
  if (!inherits(x, "STE_internal")){
    stop("Argument 'x' must be an object of class \"STE_internal\".")
  }

  no_S <- length(x$source_names)
  no_EM <- length(x$subgroup_names)

  source_num <- matrix(x$source_names, nrow = 1)
  append_row <- matrix(NA, ncol = no_S, nrow = no_EM - 1)
  source_lab <- c(rbind(source_num, append_row))

  df_dif <- x$df_dif
  df_dif$Source <- source_lab

  cat('SUBGROUP TREATMENT EFFECT ESTIMATES IN INTERNAL POPULATIONS\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(df_dif, digits = digits, ATE = FALSE, internal = TRUE)
}

#' @rdname print.STE_internal
#' @export
print.ATE_internal <- function(x, digits = 4, ...){
  if (!inherits(x, "ATE_internal")){
    stop("Argument 'x' must be an object of class \"ATE_internal\".")
  }

  cat('AVERAGE TREATMENT EFFECT ESTIMATES IN INTERNAL POPULATIONS\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(x$df_dif, digits = digits, ATE = TRUE, internal = TRUE)
}


#' @rdname print.STE_internal
#' @export
print.STE_external <- function(x, digits = 4, ...){
  if (!inherits(x, "STE_external")){
    stop("Argument 'x' must be an object of class \"STE_external\".")
  }

  cat('SUBGROUP TREATMENT EFFECT ESTIMATES IN AN EXTERNAL POPULATION\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(x$df_dif, digits = digits, ATE = FALSE, internal = FALSE)
}

#' @rdname print.STE_internal
#' @export
print.ATE_external <- function(x, digits = 4, ...){
  if (!inherits(x, "ATE_external")){
    stop("Argument 'x' must be an object of class \"ATE_external\".")
  }

  cat('AVERAGE TREATMENT EFFECT ESTIMATES IN AN EXTERNAL POPULATION\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(x$df_dif, digits = digits, ATE = TRUE, internal = FALSE)
}

#' Summary method for objects of class "ATE_internal", "ATE_external", "STE_internal", or "STE_external"
#'
#' Summary method for objects of class "ATE_internal", "ATE_external", "STE_internal", or "STE_external"
#'
#' @param object Object of class "ATE_internal", "ATE_external", "STE_internal", or "STE_external".
#' @param digits Integer specifying the number of decimal places to display.
#' @param ... Other arguments.
#' @return No value is returned.
#' @seealso \code{\link{ATE_internal}}, \code{\link{ATE_external}}, \code{\link{STE_internal}}, \code{\link{STE_external}}
#'
#' @examples
#' \donttest{
#' si <- STE_internal(
#'   X = dat_multisource[, 2:10],
#'   Y = dat_multisource$Y,
#'   EM = dat_multisource$EM,
#'   S = dat_multisource$S,
#'   A = dat_multisource$A,
#'   cross_fitting = FALSE,
#'   source_model = "MN.nnet",
#'   source_model_args = list(),
#'   treatment_model_type = "separate",
#'   treatment_model_args = list(
#'     family = binomial(),
#'     SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
#'     cvControl = list(V = 5L)
#'   ),
#'   outcome_model_args = list(
#'     family = gaussian(),
#'     SL.library = c("SL.glmnet", "SL.nnet", "SL.glm"),
#'     cvControl = list(V = 5L)
#'   )
#' )
#' summary(si)
#' }
#'
#' @export

summary.STE_internal <- function(object, digits = 4, ...){
  if (!inherits(object, "STE_internal")){
    stop("Argument 'object' must be an object of class \"STE_internal\".")
  }

  no_S <- length(object$source_names)
  no_EM <- length(object$subgroup_names)

  source_num <- matrix(object$source_names, nrow = 1)
  append_row <- matrix(NA, ncol = no_S, nrow = no_EM - 1)
  source_lab <- c(rbind(source_num, append_row))

  df_dif <- object$df_dif; df_dif$Source <- source_lab
  df_A0 <- object$df_A0; df_A0$Source <- source_lab
  df_A1 <- object$df_A1; df_A1$Source <- source_lab

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

  cat('\n\nSuperLearner libraries used:\n')
  cat("----------------------------\n")
  cat(paste0('Outcome model: ', paste0(object$outcome_model_args$SL.library, collapse = ', '), '\n'))
  cat(paste0('Treatment model: ', paste0(object$treatment_model_args$SL.library, collapse = ', '), '\n'))
  cat(paste0('Source model: NA (model fit via ', object$source_model, ')\n'))
}


#' @rdname summary.STE_internal
#' @export
summary.STE_external <- function(object, digits = 4, ...){
  if (!inherits(object, "STE_external")){
    stop("Argument 'object' must be an object of class \"STE_external\".")
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

  cat('\n\nSuperLearner libraries used:\n')
  cat("----------------------------\n")
  cat(paste0('Outcome model: ', paste0(object$outcome_model_args$SL.library, collapse = ', '), '\n'))
  cat(paste0('Treatment model: ', paste0(object$treatment_model_args$SL.library, collapse = ', '), '\n'))
  cat(paste0('Source model: NA (model fit via ', object$source_model, ')\n'))
  cat(paste0('External model: ', paste0(object$external_model_args$SL.library, collapse = ', '), '\n'))
}

#' @rdname summary.STE_internal
#' @export
summary.ATE_external <- function(object, digits = 4, ...){
  if (!inherits(object, "ATE_external")){
    stop("Argument 'object' must be an object of class \"ATE_external\".")
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

  cat('\n\nSuperLearner libraries used:\n')
  cat("----------------------------\n")
  cat(paste0('Outcome model: ', paste0(object$outcome_model_args$SL.library, collapse = ', '), '\n'))
  cat(paste0('Treatment model: ', paste0(object$treatment_model_args$SL.library, collapse = ', '), '\n'))
  cat(paste0('Source model: NA (model fit via ', object$source_model, ')\n'))
  cat(paste0('External model: ', paste0(object$external_model_args$SL.library, collapse = ', '), '\n'))
}

#' @rdname summary.STE_internal
#' @export
summary.ATE_internal <- function(object, digits = 4, ...){
  if (!inherits(object, "ATE_internal")){
    stop("Argument 'object' must be an object of class \"ATE_internal\".")
  }

  cat('AVERAGE TREATMENT EFFECT ESTIMATES IN INTERNAL POPULATIONS\n\n')
  cat('Treatment effect (mean difference) estimates:\n')
  cat("---------------------------------------------\n")
  my_print(object$df_dif, digits = digits, ATE = TRUE, internal = TRUE)

  cat('\n\nPotential outcome mean estimates under A = 0:\n')
  cat("---------------------------------------------\n")
  my_print(object$df_A0, digits = digits, ATE = TRUE, internal = TRUE)

  cat('\n\nPotential outcome mean estimates under A = 1:\n')
  cat("---------------------------------------------\n")
  my_print(object$df_A1, digits = digits, ATE = TRUE, internal = TRUE)

  cat('\n\nSuperLearner libraries used:\n')
  cat("----------------------------\n")
  cat(paste0('Outcome model: ', paste0(object$outcome_model_args$SL.library, collapse = ', '), '\n'))
  cat(paste0('Treatment model: ', paste0(object$treatment_model_args$SL.library, collapse = ', '), '\n'))
  cat(paste0('Source model: NA (model fit via ', object$source_model, ')\n'))
}



my_print <- function(df, digits, ATE, internal, ...){
  my_fun <- function(x){
    sprintf(paste0("%.", digits, "f"), round(x, digits))
  }
  mycols <- which(colnames(df) %in% c('Source', 'Subgroup'))

  if (ATE & !internal){
    df_num <- df
  } else if (ATE & internal){
    df_num <- df[, -mycols]
    df_char <- data.frame(Source = df[, mycols])
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
