glasso_gradient <- function(X_train, Y_train, X_tune = NULL, Y_tune = NULL, Xnorm,
                            lambda1_multiplier = 1, lambda, bannot,
                            no_tune = FALSE){
  if (is.null(X_tune) & is.null(Y_tune)){
    no_tune = TRUE
    X_tune = X_train
    Y_tune = Y_train
  }
  XY = sapply(1:length(Y_train), function(i){
    crossprod(X_train[[i]], Y_train[[i]]) / nrow(X_train[[i]])
  })
  XX_train = lapply(X_train, function(x){crossprod(x) / nrow(x)})
  N_train = sapply(X_train, nrow)
  mse = rep(0, nrow(lambda))

  pre_l2_est = matrix(0, ncol(X_train[[1]]), length(Y_train))
  betas = array(0, dim = c(ncol(X_train[[1]]), length(Y_train), nrow(lambda)))
  for (i in 1:nrow(lambda)){
    if (no_tune){
      res = glasso_no_early_stopping(X = X_train, Y = Y_train,
                   XX = XX_train, XY, Xnorm = Xnorm,
                   lambda1 = as.numeric(lambda[i, 1]) / N_train * lambda1_multiplier,
                   lambda2 = as.numeric(lambda[i, 2]),
                   alpha1 = as.numeric(lambda[i, 3]),
                   alpha2 = as.numeric(lambda[i, 4]),
                   theta = pre_l2_est, bannot = bannot)
      pre_l2_est = res$est
      betas[, , i] = res$est

      res$avg_tune_err = rep(0, length(X_tune))
      for (ct in 1:length(X_tune)){
        y2 = as.numeric(X_tune[[ct]] %*% res$est[, ct])
        res$avg_tune_err[ct] = mean((y2 - Y_tune[[ct]])^2)
      }
      res$avg_tune_err = mean(res$avg_tune_err)

      mse[i] = c(res$avg_tune_err)
      if (lambda[i, 1] == max(lambda[, 1])){
        pre_l2_est = matrix(0, ncol(X_train[[1]]), length(Y_train))
      }

    }else{
      res = glasso(X = X_train, Y = Y_train, X1 = X_tune, Y1 = Y_tune,
                   XX = XX_train, XY, Xnorm = Xnorm,
                   lambda1 = as.numeric(lambda[i, 1]) / N_train * lambda1_multiplier,
                   lambda2 = as.numeric(lambda[i, 2]),
                   alpha1 = as.numeric(lambda[i, 3]),
                   alpha2 = as.numeric(lambda[i, 4]),
                   theta = pre_l2_est, bannot = bannot)
      pre_l2_est = res$est
      betas[, , i] = res$est

      mse[i] = c(res$avg_tune_err)
      if (lambda[i, 1] == max(lambda[, 1])){
        pre_l2_est = matrix(0, ncol(X_train[[1]]), length(Y_train))
      }
    }
  }
  return(list(mse = mse, betas = betas))
}



#' @title cv_glasso
#'
#' @description Perform k-fold cross-validation to train a
#' @author Chen Lin, Hongyu Zhao
#' @details TBD
#' @references TBD
#' @param X_all C-list of N_c * M matrix of genotypes for training.
#' @param Y_all C-list of N_c vector of molecular phenotypes for training.
#' @param X_test C-list of N_c * M matrix of genotypes for testing.
#' @param Y_test C-list of N_c vector of molecular phenotypes for testing.
#' @param foldid N vector of fold IDs for k-fold cross-validation. The names of the vector are the individual IDs.
#' @param Xnorm M vector of the genotypic variances.
#' @param bannot M vector of the genotypic annotations. The values should be within: 1: Active; 0: Nuetral; -1: Inactive.
#' @param nlambda the number of tuning parameters lambda.
#' @param nalpha the number of tuning parameters alpha.
#' @param ... additional arguments.
#' @return A \code{"UTMOST_fa"} object with elements for the training model and testing results:
#' @importFrom dplyr filter
#' @export
cv_glasso <- function(X_all, Y_all, X_test, Y_test, foldid, Xnorm,
                      bannot, nlambda = 10, nalpha = 4, ...){
  if (min(bannot) == 1){
    bannot = bannot - 2
  } else if (min(bannot) == 0){
    bannot = bannot - 1
  }
  current_ct_list = names(Y_all)

  lambda = est_glasso_lambda(X_all, Y_all, Xnorm, bannot, nlambda1 = nlambda, nlambda2 = nlambda, nalpha = nalpha)

  unique_foldid = unique(foldid)
  nfold = length(unique_foldid)
  new_foldid = match(foldid, unique_foldid)
  names(new_foldid) = names(foldid)
  foldid = new_foldid

  mse = list()
  N = sapply(X_all, nrow)

  for (f2 in seq(nfold)){
    X_train = Y_train = X_tune = Y_tune = list()
    for (ct in current_ct_list){
      Y_train[[ct]] = Y_all[[ct]][foldid[names(Y_all[[ct]])] != f2]
      X_train[[ct]] = X_all[[ct]][foldid[names(Y_all[[ct]])] != f2, ]
      Y_tune[[ct]] = Y_all[[ct]][foldid[names(Y_all[[ct]])] == f2]
      X_tune[[ct]] = X_all[[ct]][foldid[names(Y_all[[ct]])] == f2, ]
    }

    mse[[f2]] = glasso_gradient(X_train, Y_train, X_tune, Y_tune, Xnorm,
                                lambda1_multiplier = (nfold - 1) / nfold, lambda, bannot)
  }

  lambda_min_cv = which.min(Reduce('+', lapply(mse, '[[', "mse")))
  lambda_min_cv_nfold = sapply(lapply(mse, '[[', "mse"), which.min)

  # trad re-train
  mse1 = glasso_gradient(X_all, Y_all, X_test, Y_test, Xnorm, lambda1_multiplier = 1,
                         lambda[c(lambda_min_cv, 1), ], bannot, no_tune = TRUE)
  r2p_retrain = r2p(mse1$betas[, , 1], X_test, Y_test)

  # ensemble cv
  s11 = lapply(1:nfold, function(u) mse[[u]]$betas[, , lambda_min_cv])
  beta_ensemble = Reduce("+", s11) / length(s11)
  r2p_ensemble = r2p(beta_ensemble, X_test, Y_test)

  return(list(beta_retrain = mse1$betas[, , 1],
              beta_ensemble = beta_ensemble,
              lambda = lambda,
              lambda.min = lambda_min_cv,
              r2_retrain = r2p_retrain$r2_test, p_retrain = r2p_retrain$p_test,
              r2_ensemble = r2p_ensemble$r2_test, p_ensemble = r2p_ensemble$p_test)
  )
}
