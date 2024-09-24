#' Format data for the two-stage prediction method
#'
#' Convert user input data into a data format suitable for the function sp_para().
#'
#' @param data user input data (containing instrumental variable, exposure, outcome, causal function (f), covariate (optional)).
#' @param iv_col colnames corresponding to instrumental variable.
#' @param covariate_col colnames corresponding to covariates.
#' @param exposure_col colname corresponding to exposure.
#' @param outcome_col colname corresponding to outcome.
#' @param f predefined causal function.

#' @return a data frame which is suitable for the function sp_para().
#' @export
format_sp = function(data, iv_col, covariate_col = NULL, exposure_col, outcome_col, f)

  {
  Z = as.data.frame(data[,which(colnames(data) %in% iv_col)])
  nz = ncol(Z)
  colnames(Z) = paste0('z_',1:nz)

  X = as.data.frame(data[,which(colnames(data) == exposure_col)])
  colnames(X) = 'x'

  Y = as.data.frame(data[,which(colnames(data) == outcome_col)])
  colnames(Y) = 'y'

  nf = length(f)
  for (i in 1:nf) {
    f_i = f[[i]](X)
    if (i==1) { f_d = f_i } else { f_d = cbind(f_d, f_i) }
  }
  colnames(f_d) = paste0('f_',1:nf)

  C = as.data.frame(data[,which(colnames(data) %in% covariate_col)])
  nc = ncol(C)
  if (nc != 0) {
    colnames(C) = paste0('c_',1:nc)
    d = cbind(Z,C,X,Y,f_d)
  } else {
    d = cbind(Z,X,Y,f_d)
  }

  return(d)

}

#' Conduct the two-stage prediction estimation
#'
#' conduct the two-stage prediction estimation.
#'
#' @param d output of the function format_sp()

#' @return the coefficient and corresponding standard error of every causal function f_j(X).
#' @return the F-value and corresponding P-value for hypothesis testing H_0: f = 0 (the causal effect of the exposure on the outcome is zero).
#' @export
sp_para = function(d) {

  Z = d[,which(grepl('z_', colnames(d)))]
  C = d[,which(grepl('c_', colnames(d)))]

  # stage I
  f = d[,which(grepl('f', colnames(d)))]
  nf = ncol(f)
  for (i in 1:nf) {
    if (ncol(C) != 0) {
      data = cbind(f[,i],Z,C)
    } else {
      data = cbind(f[,i],Z)
    }
    colnames(data)[1] = 'f'
    p = lm(f ~ . , data)$"fitted.values"
    p = as.data.frame(p, ncol=1)
    colnames(p)[1] = paste0('fhat_',i)
    if (i==1) { fhat = p } else { fhat = cbind(fhat, p) }
  }

  # stage II
  Y = d[,which(grepl('y', colnames(d)))]
  Y = as.data.frame(Y,ncol=1)
  colnames(Y) = 'y'

  if (ncol(C) == 0) {
    data2 = cbind(Y,fhat)
  } else {
    data2 = cbind(Y,fhat,C)
  }
  c2  = lm(y ~ ., data2)

  c2_coe = c2$coefficients
  fhat_coe = c2_coe[which(grepl('fhat_', names(c2_coe)))]

  W = model.matrix(c2)
  if (ncol(C) == 0) {
    W_real = cbind(rep(1,nrow(W)),f)
  } else {
    W_real = cbind(rep(1,nrow(W)),f,C)
  }
  delta2hat = y - as.matrix(W_real,ncol = ncol(W_real)) %*% as.matrix(c2_coe,ncol=1)
  Sigma = solve(t(W)%*%W) * as.numeric(var(delta2hat))
  Sigma_f = Sigma[2:(nf+1),2:(nf+1)]

  se_f = sqrt(diag(Sigma_f))
  names(se_f) = paste0(names(se_f),'_se')

  p_f = 2*pnorm(-abs(fhat_coe/se_f))
  names(p_f) = paste0(names(p_f),'_p')

  f_value = t(as.matrix(fhat_coe,ncol=1)) %*% solve(Sigma_f) %*% as.matrix(fhat_coe,ncol=1) / nf
  p_value = 1 - stats::pf(f_value,nf,nrow(W)-ncol(W))

  res = c(fhat_coe,se_f,p_f,f_value,p_value)
  names(res)[(length(res)-1):length(res)] = c('f_value','p_value')

  return(res)

}

