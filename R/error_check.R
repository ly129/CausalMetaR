error_check <- function(X, X_external, Y, S, A, external, ATE){

  # Checking lengths of all data inputs
  n <- nrow(X)
  if (length(Y) != n | length(S) != n | length(A) != n){
    stop(paste0('All data inputs for the sources must have the same length:\n',
                'X has length ', n, '.\n',
                'Y has length ', length(Y), '.\n',
                'S has length ', length(S), '.\n',
                'A has length ', length(A), '.'))
  }

  # Checking treatment data
  if (!is.numeric(A)){
    stop('A must be a numeric vector')
  }
  if (!all(A %in% c(0, 1))){
    stop('A must only take values 0 or 1')
  }

  # Checking source data
  if (!is.numeric(S)){
    stop('S must be a numeric vector')
  }
  unique_S <- sort(unique(S))
  no_S <- length(unique_S)
  if (!identical(seq(no_S), unique_S)) {
    stop(paste("Source", setdiff(seq(max(S)), unique_S), "is missing. "))
  }

  # Checking covariate data
  if (!ATE){
    if (!is.factor(X[, 1])){
      stop(paste0('The first column of X must be a factor variable (corresponding to the categorical potential effect modifier).'))
    }
  }
  if (external){
    if (ncol(X) != ncol(X_external)){
      stop(paste0('X and X_external must have the same number of columns.'))
    }
    if (!ATE){
      if (!is.factor(X_external[, 1])){
        stop(paste0('The first column of X_external must be a factor variable (corresponding to the categorical potential effect modifier).'))
      }
    }
  }
}
