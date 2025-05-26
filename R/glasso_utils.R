glasso_obj <- function(theta, bannot, Y, X, lambda1, lambda2){
  P = ncol(theta)
  theta_L2 = sum(sqrt(rowSums(theta^2)) * lambda2 * bannot)

  new_objV1_train = new_objV1_obj = rep(0, P)
  for (t in 1:P){
    new_objV1_train[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
    new_objV1_obj[t] =  new_objV1_train[t] + lambda1[t] * sum(abs(theta[, t])) + theta_L2 / P
  }
  return(cbind(new_objV1_train, new_objV1_obj))
}

#' @importFrom dplyr filter %>%
est_glasso_lambda <- function(X_all, Y_all, Xnorm, bannot,
                              nlambda1 = 10, nlambda2 = 10, nalpha = 4){
  XY = sapply(1:length(Y_all), function(i){
    crossprod(X_all[[i]], Y_all[[i]]) / nrow(X_all[[i]])
  })
  N = sapply(X_all, nrow)

  lambda_ct = sapply(1:length(X_all), function(ct)
    crossprod(X_all[[ct]], Y_all[[ct]]) %>% abs %>% max)
  lambda_max = max(lambda_ct)
  lambda_min = min(lambda_ct) * 0
  lambda1 = seq(lambda_max, lambda_min, length.out = nlambda2 * 5)

  dense_seq = (nlambda2 * 5 + 1) - ((seq(2, nlambda2 + 1) / (nlambda2+1))^2 * (nlambda2 * 5 - 4)) %>% round %>% unique
  lambda2 = sapply(dense_seq, function(i) {
    beta_lasso = pmax(sweep(abs(XY), 2, lambda1[i] / N, "-"), 0) / Xnorm
    theta_L2 = beta_lasso^2 %>% rowSums %>% sqrt
    max(theta_L2)
  })
  lambda1_max = lambda1[dense_seq]

  lambda_list = list()
  for (i in 1:nlambda2){
    lambda1 = seq(lambda1_max[i], lambda1_max[i] * 0.05, length.out = nlambda1 + 1) %>% exp
    lambda1 = lambda1[-1]

    alpha_grid = expand.grid(lambda1, lambda2[i],
                             ((nalpha):1) / (nalpha),
                             ((nalpha):1) / (nalpha))
    alpha_grid <- filter(alpha_grid, Var3 <= Var4)
    lambda_list[[i]] = alpha_grid
  }
  lambda = do.call(rbind, lambda_list)
  return(lambda)
}

##############################################################
## group lasso on data with missing covariates, use validation data for stopping criterion
glasso <- function(X, Y, X1, Y1, XX, XY, Xnorm,
                   lambda1, lambda2, alpha1, alpha2, theta, bannot,
                   maxiter = 50, eps = 1e-3, verbose = FALSE){
  M = nrow(XY)
  P = length(X)
  NN = unlist(lapply(X, nrow))
  bannot = ifelse(bannot == 1, alpha1, 0) + ifelse(bannot == 0, alpha2, 0) + ifelse(bannot == -1, 1, 0)

  old_objV1 = glasso_obj(theta, bannot, Y, X, lambda1, lambda2)[, 2]
  old_objV2 = rep(0,P)
  for(t in 1:P){
    old_objV2[t] = mean((Y1[[t]]-X1[[t]]%*%theta[,t])^2)
  }
  beta_j_lasso = rep(0, P)

  # if(!is.loaded("wrapper2")){
  #   dyn.load("/gpfs/gibbs/pi/zhao/cl2384/CONDA/conda_envs/CTIMP/optima.so") # change this to the abs path to optim.so
  # }
  for(i in 1:maxiter){
    res = .Call("wrapper2", XX, XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm, bannot, PACKAGE = "UTMOSTfa")

    new_objV1 = glasso_obj(theta, bannot, Y, X, lambda1, lambda2)[, 2]
    new_objV2 = rep(0,P)
    for(t in 1:P){
      new_objV2[t] = mean((Y1[[t]]-X1[[t]]%*%theta[,t])^2)
    }

    if (mean(new_objV2) >= mean(old_objV2) | mean(new_objV1) >= mean(old_objV1) | max(abs(new_objV1 - old_objV1)) < eps){
      break
    }

    old_objV2 = new_objV2
    old_objV1 = new_objV1

  }
  list(est = theta, avg_tune_err = mean(new_objV2), tune_err=new_objV2)
  #avg_obj_err = mean(new_objV1), avg_train_err = mean(glasso_obj(theta, bannot, Y, X, lambda1, lambda2)[, 1]))
}

## simpler version of glasso, train model until converges
glasso_no_early_stopping <- function(X, Y, XX, XY, Xnorm,
                                     lambda1, lambda2, alpha1, alpha2, theta, bannot,
                                     maxiter = 50, eps = 1e-3, verbose = FALSE){
  M = nrow(XY)
  P = length(X)
  NN = unlist(lapply(X, nrow))
  bannot = ifelse(bannot == 1, alpha1, 0) + ifelse(bannot == 0, alpha2, 0) + ifelse(bannot == -1, 1, 0)

  old_objV1 = glasso_obj(theta, bannot, Y, X, lambda1, lambda2)[, 2]

  beta_j_lasso = rep(0, P)
  # if(!is.loaded("wrapper2")){
  #   dyn.load("/gpfs/gibbs/pi/zhao/cl2384/CONDA/conda_envs/CTIMP/optima.so")
  # }
  for(i in 1:maxiter){
    res = .Call("wrapper2", XX, XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm, bannot, PACKAGE = "UTMOSTfa")

    new_objV1 = glasso_obj(theta, bannot, Y, X, lambda1, lambda2)[, 2]
    if(max(abs(new_objV1-old_objV1)) < eps | mean(new_objV1) > mean(old_objV1)){
      break
    }
    old_objV1 = new_objV1
  }
  list(est = theta, avg_obj_err = mean(new_objV1),
       avg_train_err = mean(glasso_obj(theta, bannot, Y, X, lambda1, lambda2)[, 1]), obj_err = new_objV1)
}
