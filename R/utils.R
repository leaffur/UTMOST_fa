ACAT_test <- function(r2, p){
  T_ACAT = mean(tan((0.5 - p) * pi))
  p_ACAT = 0.5 - atan(T_ACAT) / pi
  return(c(mean(r2), p_ACAT))
}

#' @importFrom stats cor.test var
r2p <- function(est, X_test, Y_test){
  current_ct_list = names(Y_test)
  r2_test = rep(0, length(current_ct_list))
  p_test = rep(1, length(current_ct_list))
  names(r2_test) = names(p_test) = current_ct_list

  if (is.null(dim(est))){
    est = matrix(est, ncol = 1)
  }
  for (ct in 1:length(current_ct_list)){
    y2 = as.numeric(X_test[[ct]] %*% est[, ct])
    if (var(y2) > 0 & var(Y_test[[ct]]) > 0){
      c1 = cor.test(y2, Y_test[[ct]])
      r2_test[ct] = c1$estimate^2
      p_test[ct] = c1$p.value
    }
  }
  return(list(r2_test = r2_test, p_test = p_test))
}

#' @importFrom dplyr %>%
est_lasso_lambda1 <- function(X_all, Y_all, nlambda = 50){
  current_ct_list = names(Y_all)
  lambda = list()

  for (ct in current_ct_list){
    lambda_max = crossprod(X_all[[ct]], Y_all[[ct]]) %>% abs %>% max / length(Y_all[[ct]])
    lambda.ratio = 0.1^(2 / 100 * nlambda)
    lambda[[ct]] = seq(log(lambda_max), log(lambda_max * lambda.ratio), length.out = nlambda) %>% exp %>% round(digits = 8)
  }
  return(lambda)
}

#' @importFrom glmnet cv.glmnet
#' @importFrom dplyr %>%
#' @importFrom stats predict coef
cv_glmnet <- function(X_all, Y_all, X_test, Y_test, foldid, nlambda = 100){
  t1 = Sys.time()
  lambda = est_lasso_lambda1(X_all, Y_all, nlambda)

  current_ct_list = names(Y_all)
  r2_test_lasso = lambda.min = rep(0, length(current_ct_list))
  p_test_lasso = rep(1, length(current_ct_list))
  names(r2_test_lasso) = names(lambda.min) = names(p_test_lasso) = current_ct_list

  unique_foldid = unique(foldid)

  beta = matrix(0, 1 + ncol(X_all[[1]]), length(current_ct_list))
  colnames(beta) = current_ct_list

  cvm = rep(0, nlambda)
  for (ct in current_ct_list){
    foldid_new = foldid[match(names(Y_all[[ct]]), names(foldid))]
    foldid_new <- match(foldid_new, unique_foldid)
    fit2 = cv.glmnet(X_all[[ct]], Y_all[[ct]], alpha = 0.5,
                     standardize = FALSE,
                     foldid = foldid_new,
                     lambda = lambda[[ct]] * 2, parallel = FALSE)
    cvm = cvm + fit2$cvm
    y2 = predict(fit2, newx = X_test[[ct]], s = "lambda.min")
    lambda.min[ct] = fit2$lambda.min
    beta[, ct] = coef(fit2, s = "lambda.min") %>% as.numeric # c(fit2$glmnet.fit$a0[fit2$index[1]], fit2$glmnet.fit$beta[, fit2$index[1]])
    if (var(y2) > 0 & var(Y_test[[ct]]) > 0){
      c1 = cor.test(y2, Y_test[[ct]])
      r2_test_lasso[ct] = c1$estimate^2
      p_test_lasso[ct] = c1$p.value
    }
  }

  t2 = Sys.time()
  return(list(beta = beta,
              cvm = cvm / length(current_ct_list),
              lambda.min = lambda.min,
              r2_test = r2_test_lasso,
              p_test = p_test_lasso))
}
