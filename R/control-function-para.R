#' Format data for the control function method
#'
#' Convert user input data into a data format suitable for the function cf_para().
#'
#' @param data user input data (containing instrumental variable, exposure, outcome, f function, covariate(optional), g function(optional)).
#' @param iv_col colnames corresponding to instrumental variables.
#' @param covariate_col colnames corresponding to covariates (optional).
#' @param exposure_col colname corresponding to exposure.
#' @param outcome_col colname corresponding to outcome.
#' @param f predefined causal function of exposure on outcome.
#' @param g predefined causal function of covariate on outcome (optional).

#' @return a data frame which is suitable for the function cf_para().
#' @export
format_cf = function(data, iv_col, covariate_col = NULL, exposure_col, outcome_col, f, g = NULL)

  {

  Z = as.data.frame(data[,which(colnames(data) %in% iv_col)])
  nz = ncol(Z)
  colnames(Z) = paste0('z_',1:nz)

  C = as.data.frame(data[,which(colnames(data) %in% covariate_col)])
  nc = ncol(C)
  if (nc != 0) {
    colnames(C) = paste0('c_',1:nc)
  }

  X = data[,which(colnames(data) == exposure_col)]
  X = as.data.frame(X,ncol=1)
  colnames(X) = 'x'

  Y = data[,which(colnames(data) == outcome_col)]
  Y = as.data.frame(Y,ncol=1)
  colnames(Y) = 'y'

  nf = length(f)
  for (i in 1:nf) {
    f_i = f[[i]](X)
    if (i==1) { f_d = f_i } else { f_d = cbind(f_d, f_i) }
  }
  f_d = as.data.frame(f_d, ncol=nf)
  colnames(f_d) = paste0('f_',1:nf)

  if (!is.null(g)) {
    g_d = as.data.frame(g(C), ncol=1)
    names(g_d) = 'g'
  }

  if (is.null(g)) {
    d = cbind(Z,X,Y,f_d)
  } else {
    d = cbind(Z,C,X,Y,f_d,g_d)
  }

  return(d)

}

#' var_est_gaussian_h0
#'
#' A function used in function cf_para().
#' @export
var_est_gaussian_h0 = function(c1,c2) {
  W = model.matrix(c2)
  V = model.matrix(c1)
  c2_coe = c2$coefficients
  rhohat = c2_coe[which(grepl('residuals', names(c2_coe)))][[1]]
  ehat = c2$residuals
  V_beta = vcov(c1)
  v2 = vcov(c2)
  sigma = v2 + (v2/var(ehat)) %*% (rhohat^2 * t(W) %*% V %*% V_beta %*% t(V) %*% W ) %*% (v2/var(ehat))
  return(sigma)
}


#' var_est_gaussian_h1
#'
#' A function used in function cf_para().
#' @export
var_est_gaussian_h1 = function(c1,c2,h,h_d) {
  W = model.matrix(c2)
  V = model.matrix(c1)
  c2_coe = c2$coefficients
  ehat = c2$residuals
  evv = nrow(W)*solve(t(V)%*%V)
  delta1hat = c1$residuals
  rhohat = c2_coe[which(grepl('residuals', names(c2_coe)))][[1]]
  W_n = cbind(rhohat*W[,1:(ncol(W)-1)], rhohat*W[,ncol(W)]-ehat)
  de = matrix(rep(NA,ncol(W_n)*ncol(V)), nrow = ncol(V), ncol = ncol(W_n))
  for (i in 1:ncol(V)) {
    for (j in 1:ncol(W_n)) {
      de[i,j] = mean(h_d(delta1hat) * V[,i] * W_n[,j])
    }
  }
  A = W * ehat + matrix(V*delta1hat,ncol=ncol(V)) %*% (evv %*% de)
  Sigma = t(A)%*%A/nrow(W)
  sigma = nrow(W)*solve(t(W)%*%W) %*% Sigma %*% (nrow(W)*solve(t(W)%*%W)) / nrow(W)

  return(sigma)

}


#' var_est_binary_h0
#'
#' A function used in function cf_para().
#' @export
var_est_binary_h0 = function(c1,c2,data2) {
  W = model.matrix(c2)
  V = model.matrix(c1)
  c2_coe = c2$coefficients
  rhohat = c2_coe[which(grepl('residuals', names(c2_coe)))][[1]]
  r2 = matrix(c2$residuals,ncol=1)
  etahat = c2$linear.predictors
  miuhat = c2$fitted.values
  y = data2$y
  gz = W * as.numeric(y-miuhat)
  fiez = as.matrix(V * c1$residuals, ncol=ncol(V)) %*% solve(crossprod(V)) * nrow(data2)
  W_n = cbind(rhohat*miuhat*(1-miuhat)*W[,1:(ncol(W)-1)],
              (c1$residuals*rhohat*miuhat*(1-miuhat)-(y-miuhat)))
  de = matrix(rep(NA,ncol(W_n)*ncol(V)), nrow = ncol(V), ncol = ncol(W_n))
  for (i in 1:ncol(V)) {
    for (j in 1:ncol(W_n)) {
      de[i,j] = mean(as.data.frame(V)[,i] * W_n[,j])
    }
  }
  A = gz + fiez %*% de
  v = vcov(c2)
  sigma = v %*% crossprod(A) %*% v
  return(sigma)
}


#' var_est_binary_h1
#'
#' A function used in function cf_para().
#' @export
var_est_binary_h1 = function(c1,c2,data2,h,h_d) {
  W = model.matrix(c2)
  V = model.matrix(c1)
  c2_coe = c2$coefficients
  rhohat = c2_coe[which(grepl('residuals', names(c2_coe)))][[1]]
  r2 = matrix(c2$residuals,ncol=1)
  etahat = c2$linear.predictors
  miuhat = c2$fitted.values
  delta1hat = c1$residuals
  y = data2$y
  gz = W * as.numeric(y-miuhat)
  fiez = as.matrix(V * c1$residuals, ncol=ncol(V)) %*% solve(crossprod(V)) * nrow(data2)
  W_n = cbind(rhohat*miuhat*(1-miuhat)*W[,1:(ncol(W)-1)], (h(delta1hat)*rhohat*miuhat*(1-miuhat)-(y-miuhat)))
  de = matrix(rep(NA,ncol(W_n)*ncol(V)), nrow = ncol(V), ncol = ncol(W_n))
  for (i in 1:ncol(V)) {
    for (j in 1:ncol(W_n)) {
      de[i,j] = mean(h_d(delta1hat)[,1] * as.data.frame(V)[,i] * W_n[,j])
    }
  }
  A = gz + fiez %*% de
  v = vcov(c2)
  sigma = v %*% crossprod(A) %*% v
  return(sigma)
}


#' Control function method
#'
#' conduct the control function estimation.
#'
#' @param d output of the function format_cf()
#' @param family c('gaussian','binary'): 'gaussian' for the normal outcome variable; 'binary' for the binary outcome variable.
#' @param hp c(0,1): 0 for the scenario where the horizontal pleiotropy is absent and 1 for the scenario where the horizontal pleiotropy is present.
#' @param h the function h(delta1), the default is a linear function.
#' @param h_d the derivative of function h(delta1).

#' @return the coefficient, corresponding standard error and P-value of every causal function f_j(X).
#' @return the F-value and corresponding P-value for hypothesis testing H_0: f = 0 (the causal effect of the exposure on the outcome is zero).
#' @export
cf_para = function(d,
                   family = 'gaussian',    # c('gaussian','binary')
                   hp = 0,  # whether the horizontal pleiotropy is present; 1: present; 0: absent
                   h = NULL,
                   h_d = NULL
                   )
{

  # d: output of format_cf()

  # Stage I - Data generation
  Z = d[,which(grepl('z_', colnames(d)))]
  C = d[,which(grepl('c_', colnames(d)))]
  if (length(C) == 0) {
    data = as.data.frame(cbind(d$x,Z))
  } else {
    data = as.data.frame(cbind(d$x,Z,C))
  }
  names(data)[1] = 'x'
  # Stage I - Estimation
  c1 = lm(x ~ ., data)
  delta1hat = as.data.frame(c1$residuals, ncol=1)

  # Stage II - Data generation
  Y = d[,which(grepl('y', colnames(d)))]
  Y = as.data.frame(Y,ncol=1)
  colnames(Y) = 'y'

  f = d[,which(grepl('f_', colnames(d)))]
  f = as.data.frame(f,nrow = nrow(f))

  if (length(C) == 0) {
    data2 = cbind(Y,f)
  } else {
    data2 = cbind(Y,f,d$g)
  }

  if (is.null(h) & is.null(h_d)) {
    data2 = cbind(data2, delta1hat)
  } else {
    data2 = cbind(data2,h(delta1hat))
  }

  if (hp == 1) {
    data2 = cbind(data2,Z)
  }

  # Stage II - Estimation
  if (family == 'gaussian') {
    c2  = lm(y ~ ., data2)
    # coe
    c2_coe = c2$coefficients
    f_coe = c2_coe[which(grepl('f_', names(c2_coe)))]
    # asymptotic variance
    if (is.null(h) & is.null(h_d)) {
      sigma = spmr::var_est_gaussian_h0(c1,c2)
    } else {
      sigma = spmr::var_est_gaussian_h1(c1,c2,h,h_d)
    }
  } else if (family == 'binary') {
    c2  = glm(y ~ ., data = data2, family = "binomial")
    # coe
    c2_coe = c2$coefficients
    f_coe = c2_coe[which(grepl('f', names(c2_coe)))]
    # asymptotic variance
    if (is.null(h) & is.null(h_d)) {
      sigma = spmr::var_est_binary_h0(c1,c2,data2)
    } else {
      sigma = spmr::var_est_binary_h1(c1,c2,data2,h,h_d)
    }
  }

  Sigma_f = sigma[2:(ncol(f)+1),2:(ncol(f)+1)]
  if (ncol(f) > 1) {
    se_f = sqrt(diag(Sigma_f))
    names(se_f) = paste0(names(se_f),'_se')
  } else if (ncol(f) == 1) {
    se_f = sqrt(Sigma_f)
    names(se_f) = paste0('f','_se')
  }

  p_f = 2*pnorm(-abs(f_coe/se_f))
  names(p_f) = paste0(names(p_f),'_p')

  f_value = t(as.matrix(f_coe,ncol=1)) %*% solve(Sigma_f) %*% as.matrix(f_coe,ncol=1) / ncol(f)
  W = model.matrix(c2)
  p_value = 1 - stats::pf(f_value,ncol(f),nrow(W)-ncol(W))

  res = c(f_coe,se_f,p_f,f_value,p_value)
  names(res)[(length(res)-1):length(res)] = c('f_value','p_value')

  return(res)

}

