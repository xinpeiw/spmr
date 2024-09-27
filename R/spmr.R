#' var_cal_binomial_h
#'
#' A function used in function cf_semipara().
#' @export
var_cal_binomial_h = function(object1,object2,y,h_d_est)
  # object1: stage 1 model
  # object2: gam model
{
  W = model.matrix(object2)
  c2_coe = object2$coefficients
  rhohat = c2_coe[which(grepl("h_r", names(c2_coe)))][[1]]
  B = matrix(c2_coe,ncol=1)
  etahat = object2$linear.predictors
  miuhat = object2$fitted.values
  gz = W * as.numeric(y-miuhat)

  V = model.matrix(object1)
  delta1hat = object1$residuals
  N = length(y)
  fiez = as.matrix(V * delta1hat, ncol=ncol(V)) %*% solve(crossprod(V)) * N
  W_n = cbind(rhohat*miuhat*(1-miuhat)*W[,1:2],
              rhohat*W[,3]*miuhat*(1-miuhat)-(y-miuhat),
              rhohat*miuhat*(1-miuhat)*W[,4:12])
  de = matrix(rep(NA,ncol(W_n)*ncol(V)), nrow = ncol(V), ncol = ncol(W_n))
  for (i in 1:ncol(V)) {
    for (j in 1:ncol(W_n)) {
      de[i,j] = mean(h_d_est(delta1hat) * as.data.frame(V)[,i] * W_n[,j])
    }
  }
  A = gz + fiez %*% de
  S = (object2$smooth[[1]])$S[[1]]
  S = rbind(matrix(rep(0,ncol(W)*(ncol(W)-nrow(S))),nrow = (ncol(W)-nrow(S))),
            cbind(matrix(rep(0,(ncol(W)-nrow(S))*nrow(S)),ncol = (ncol(W)-nrow(S))),S))
  lambda = object2$sp[[1]]
  sigma = vcov(object2) %*% crossprod(A) %*% solve(solve(vcov(object2))-lambda*S)
  return(sigma)
}


#' var_cal_binomial_r
#'
#' A function used in function cf_semipara().
#' @export
var_cal_binomial_r = function(object1,object2,y)
  # object1: stage 1 model
  # object2: gam model
{
  W = model.matrix(object2)
  c2_coe = object2$coefficients
  rhohat = c2_coe[which(grepl("h_r", names(c2_coe)))][[1]]
  B = matrix(c2_coe,ncol=1)
  etahat = object2$linear.predictors
  miuhat = object2$fitted.values
  gz = W * as.numeric(y-miuhat)

  V = model.matrix(object1)
  delta1hat = object1$residuals
  N = length(y)
  fiez = as.matrix(V * delta1hat, ncol=ncol(V)) %*% solve(crossprod(V)) * N
  W_n = cbind(rhohat*miuhat*(1-miuhat)*W[,1:2],
              rhohat*delta1hat*miuhat*(1-miuhat)-(y-miuhat),
              rhohat*miuhat*(1-miuhat)*W[,4:12])
  de = matrix(rep(NA,ncol(W_n)*ncol(V)), nrow = ncol(V), ncol = ncol(W_n))
  for (i in 1:ncol(V)) {
    for (j in 1:ncol(W_n)) {
      de[i,j] = mean(as.data.frame(V)[,i] * W_n[,j])
    }
  }
  A = gz + fiez %*% de
  S = (object2$smooth[[1]])$S[[1]]
  S = rbind(matrix(rep(0,ncol(W)*(ncol(W)-nrow(S))),nrow = (ncol(W)-nrow(S))),
            cbind(matrix(rep(0,(ncol(W)-nrow(S))*nrow(S)),ncol = (ncol(W)-nrow(S))),S))
  lambda = object2$sp[[1]]
  sigma = vcov(object2) %*% crossprod(A) %*% solve(solve(vcov(object2))-lambda*S)
  return(sigma)
}


#' var_cal_gaussian_r
#'
#' A function used in function cf_semipara().
#' @export
var_cal_gaussian_r = function(object1,object2,y)
  # object1: stage 1 model
  # object2: gam model
{
  W = predict(object2,type='lpmatrix')
  c2 = lm(y ~ W-1)
  rhohat = c2$coefficients[[which(names(c2$coefficients)=='Wh_r')]]
  ehat = c2$residuals
  delta1hat = object1$residuals
  V = model.matrix(object1)
  V_beta = vcov(object1,freq = TRUE)
  S = (object2$smooth[[1]])$S[[1]]
  S = rbind(matrix(rep(0,ncol(W)*(ncol(W)-nrow(S))),nrow = (ncol(W)-nrow(S))),
            cbind(matrix(rep(0,(ncol(W)-nrow(S))*nrow(S)),ncol = (ncol(W)-nrow(S))),S))
  lambda = object2$sp[[1]]
  sigma = solve(crossprod(W)+lambda*S) * var(ehat) +
    solve(crossprod(W)+lambda*S) %*% (rhohat^2 * t(W) %*% V %*% V_beta %*% t(V) %*% W ) %*% solve(crossprod(W))

  return(sigma)
}


#' var_cal_gaussian_h
#'
#' A function used in function cf_semipara().
#' @export
var_cal_gaussian_h = function(object1,object2,y,h_d_est)
  # object1: stage 1 model
  # object2: gam model
{
  W = predict(object2,type='lpmatrix')
  N = nrow(W)
  c2 = lm(y ~ W-1)
  rhohat = c2$coefficients[[which(names(c2$coefficients)=='Wh_r')]]
  ehat = c2$residuals
  delta1hat = object1$residuals
  V = model.matrix(object1)
  V_beta = vcov(object1,freq = TRUE)
  W_n = cbind(rhohat*W[,1:2],
              rhohat*W[,3]-ehat,
              rhohat*W[,4:12])
  de = matrix(rep(NA,ncol(W_n)*ncol(V)), nrow = ncol(V), ncol = ncol(W_n))
  for (i in 1:ncol(V)) {
    for (j in 1:ncol(W_n)) {
      de[i,j] = mean(h_d_est(delta1hat) * V[,i] * W_n[,j])
    }
  }
  A = W * ehat + matrix(V*delta1hat,ncol=ncol(V)) %*% (N*solve(crossprod(V)) %*% de)
  S = (object2$smooth[[1]])$S[[1]]
  S = rbind(matrix(rep(0,ncol(W)*(ncol(W)-nrow(S))),nrow = (ncol(W)-nrow(S))),
            cbind(matrix(rep(0,(ncol(W)-nrow(S))*nrow(S)),ncol = (ncol(W)-nrow(S))),S))
  lambda = object2$sp[[1]]
  sigma = solve(crossprod(W)+lambda*S) %*% crossprod(A) %*% solve(crossprod(W))
  return(sigma)
}



#' testStat
#'
#' A function used in function cf_semipara().
#' @export
testStat <- function(p,X,V,rank=NULL,type=0,res.df= -1) {
  ## Implements Wood (2013) Biometrika 100(1), 221-228
  ## The type argument specifies the type of truncation to use.
  ## on entry `rank' should be an edf estimate
  ## 0. Default using the fractionally truncated pinv.
  ## 1. Round down to k if k<= rank < k+0.05, otherwise up.
  ## res.df is residual dof used to estimate scale. <=0 implies
  ## fixed scale.

  qrx <- qr(X,tol=0)
  R <- qr.R(qrx)
  V <- R%*%V[qrx$pivot,qrx$pivot,drop=FALSE]%*%t(R)
  V <- (V + t(V))/2
  ed <- eigen(V,symmetric=TRUE)
  ## remove possible ambiguit y from statistic...
  siv <- sign(ed$vectors[1,]);siv[siv==0] <- 1
  ed$vectors <- sweep(ed$vectors,2,siv,"*")

  k <- max(0,floor(rank))
  nu <- abs(rank - k)     ## fractional part of supplied edf
  if (type==1) { ## round up is more than .05 above lower
    if (rank > k + .05||k==0) k <- k + 1
    nu <- 0;rank <- k
  }

  if (nu>0) k1 <- k+1 else k1 <- k

  ## check that actual rank is not below supplied rank+1
  r.est <- sum(ed$values > max(ed$values)*.Machine$double.eps^.9)
  if (r.est<k1) {k1 <- k <- r.est;nu <- 0;rank <- r.est}

  ## Get the eigenvectors...
  # vec <- qr.qy(qrx,rbind(ed$vectors,matrix(0,nrow(X)-ncol(X),ncol(X))))
  vec <- ed$vectors
  if (k1<ncol(vec)) vec <- vec[,1:k1,drop=FALSE]

  ## deal with the fractional part of the pinv...
  if (nu>0&&k>0) {
    if (k>1) vec[,1:(k-1)] <- t(t(vec[,1:(k-1)])/sqrt(ed$val[1:(k-1)]))
    b12 <- .5*nu*(1-nu)
    if (b12<0) b12 <- 0
    b12 <- sqrt(b12)
    B <- matrix(c(1,b12,b12,nu),2,2)
    ev <- diag(ed$values[k:k1]^-.5,nrow=k1-k+1)
    B <- ev%*%B%*%ev
    eb <- eigen(B,symmetric=TRUE)
    rB <- eb$vectors%*%diag(sqrt(eb$values))%*%t(eb$vectors)
    vec1 <- vec
    vec1[,k:k1] <- t(rB%*%diag(c(-1,1))%*%t(vec[,k:k1]))
    vec[,k:k1] <- t(rB%*%t(vec[,k:k1]))
  } else {
    vec1 <- vec <- if (k==0) t(t(vec)*sqrt(1/ed$val[1])) else
      t(t(vec)/sqrt(ed$val[1:k]))
    if (k==1) rank <- 1
  }
  ## there is an ambiguity in the choise of test statistic, leading to slight
  ## differences in the p-value computation depending on which of 2 alternatives
  ## is arbitrarily selected. Following allows both to be computed and p-values
  ## averaged (can't average test stat as dist then unknown)
  d <- t(vec)%*%(R%*%p)
  d <- sum(d^2)
  d1 <- t(vec1)%*%(R%*%p)
  d1 <- sum(d1^2)
  ##d <- d1 ## uncomment to avoid averaging

  rank1 <- rank ## rank for lower tail pval computation below

  ## note that for <1 edf then d is not weighted by EDF, and instead is
  ## simply refered to a chi-squared 1

  if (nu>0) { ## mixture of chi^2 ref dist
    if (k1==1) rank1 <- val <- 1 else {
      val <- rep(1,k1) ##ed$val[1:k1]
      rp <- nu+1
      val[k] <- (rp + sqrt(rp*(2-rp)))/2
      val[k1] <- (rp - val[k])
    }
    library(mgcv)
    if (res.df <= 0) pval <- (psum.chisq(d,val)+psum.chisq(d1,val))/2 else {  ## (liu2(d,val) + liu2(d1,val))/2 else
      k0 <- max(1,round(res.df))
      pval <- (psum.chisq(0,c(val,-d/k0),df=c(rep(1,length(val)),k0)) + psum.chisq(0,c(val,-d1/k0),df=c(rep(1,length(val)),k0)) )/2
      #pval <- (simf(d,val,res.df) + simf(d1,val,res.df))/2
    }
  } else { pval <- 2 }
  ## integer case still needs computing,
  ## OLD: also liu/pearson approx only good in
  ## upper tail. In lower tail, 2 moment approximation is better (Can check this
  ## by simply plotting the whole interesting range as a contour plot!)
  ##if (pval > .5)
  if (pval > 1) {
    if (res.df <= 0) pval <- (pchisq(d,df=rank1,lower.tail=FALSE)+pchisq(d1,df=rank1,lower.tail=FALSE))/2 else
      pval <- (pf(d/rank1,rank1,res.df,lower.tail=FALSE)+pf(d1/rank1,rank1,res.df,lower.tail=FALSE))/2
  }
  list(stat=d,pval=min(1,pval),rank=rank)
} ## end of testStat



#' p_cal
#'
#' A function used in function cf_semipara().
#' @export
p_cal = function(object,V,i)
  # object: gam; V: variance; i: No.smoothing term
{
  start <- object$smooth[[i]]$first.para;stop <- object$smooth[[i]]$last.para
  V <- V[start:stop,start:stop,drop=FALSE]
  p <- object$coefficients[start:stop]  # params for smooth
  edf1i <- edfi <- sum(object$edf[start:stop]) # edf for this smooth
  ## extract alternative edf estimate for this smooth, if possible...
  if (!is.null(object$edf1)) edf1i <-  sum(object$edf1[start:stop])
  if (is.null(object$R)) { ## Factor from QR decomp of sqrt(W)X
    warning("p-values for any terms that can be penalized to zero will be unreliable: refit model to fix this.")
    useR <- FALSE
  } else useR <- TRUE

  if (useR)  X <- object$R else {
    X <- model.matrix(object)
  }
  X <- X[!is.na(rowSums(X)),]
  Xt <- X[,start:stop,drop=FALSE]
  fx <- if (inherits(object$smooth[[i]],"tensor.smooth")&&
            !is.null(object$smooth[[i]]$fx)) all(object$smooth[[i]]$fx) else object$smooth[[i]]$fixed
  est.disp <- object$scale.estimated
  residual.df = length(object$y)-sum(object$edf)
  if (!fx&&object$smooth[[i]]$null.space.dim==0&&!is.null(object$R)) { ## random effect or fully penalized term
    res <- if (re.test) reTest(object,i) else NULL
  } else { ## Inverted Nychka interval statistics

    if (est.disp) rdf <- residual.df else rdf <- -1
    res <- testStat(p,Xt,V,min(ncol(Xt),edf1i),type=0,res.df = rdf)
  }

  return(res$pval)
}



#' Format data for the semi-parametric control function method
#'
#' Convert user input data into a data format suitable for the function cf_semipara().
#'
#' @param data user input data (containing instrumental variable, exposure, outcome, covariate(optional)).
#' @param iv_col colnames corresponding to instrumental variables.
#' @param covariate_col colnames corresponding to covariates (optional).
#' @param exposure_col colname corresponding to exposure.
#' @param outcome_col colname corresponding to outcome.


#' @return a data frame which is suitable for the function cf_semipara().
#' @export
format_cf_semipara = function(data, iv_col, covariate_col = NULL, exposure_col, outcome_col)

{

  # data: input data ( containing iv, exposure, outcome, f function, covariate(optional) )
  # iv_col: colnames corresponding to IV
  # covariate_col: colnames corresponding to covariates
  # exposure_col: colnames corresponding to exposure
  # outcome_col: colnames corresponding to outcome

  Z = as.data.frame(data[,which(colnames(data) %in% iv_col)])
  nz = ncol(Z)
  colnames(Z) = paste0('z_',1:nz)

  X = data[,which(colnames(data) == exposure_col)]
  X = as.data.frame(X,ncol=1)
  colnames(X) = 'x'

  Y = data[,which(colnames(data) == outcome_col)]
  Y = as.data.frame(Y,ncol=1)
  colnames(Y) = 'y'

  C = as.data.frame(data[,which(colnames(data) %in% covariate_col)])
  nc = ncol(C)
  if (nc != 0) {
    colnames(C) = paste0('c_',1:nc)
    d = cbind(Z,C,X,Y)
  } else {
    d = cbind(Z,X,Y)
  }

  return(d)

}


#' Semi-parametric control function method
#'
#' conduct the semi-parametric control function estimation.
#'
#' @param d output of the function format_cf_semipara().
#' @param family c('gaussian','binary'): 'gaussian' for the normal outcome variable; 'binary' for the binary outcome variable.
#' @param g_est predefined g function (the causal effect of covariate on the outcome), and the default is linear function.
#' @param action 0: the relationship between delta2 and delta1 is linear.
#' @param h_est the relationship between delta2 and delta1. if action!=0, h_est should be defined.
#' @param h_d_est the derivative of function h_est(). If action!=0, h_d_est should be defined.
#'
#'
#' @return the P-value for hypothesis testing H_0: f = 0 (the causal effect of the exposure on the outcome is zero).
#' @export
cf_semipara = function(d,
                       family = NULL,    # c('gaussian','binary')
                       g_est = function(d) {return(d)},  # g should be defined
                       action = 0,  # 0: the relationship between delta2 and delta1 is linear
                       h_est = NULL, # if action!=0, h_est should be defined
                       h_d_est = NULL # if action!=0, h_d_est should be defined

) {
  # d: output of format_cf_semipara ()

  Z = d[,which(grepl('z_', colnames(d)))]
  C = d[,which(grepl('c_', colnames(d)))]

  # stage I -- Data generation
  if (length(C)==0) {
    data = as.data.frame(cbind(d$x,Z))
  } else {
    data = as.data.frame(cbind(d$x,Z,C))
  }
  names(data)[1] = 'x'
  # Stage I -- Estimation
  c1 = lm(x ~ ., data)
  r = as.data.frame(c1$residuals, ncol=1)

  # stage II -- Data generation
  Y = d[,which(grepl('y', colnames(d)))]
  Y = as.data.frame(Y,ncol=1)
  colnames(Y) = 'y'
  X = d[,which(grepl('x', colnames(d)))]
  X = as.data.frame(X,ncol=1)
  colnames(X) = 'x'

  if (action == 0) {
    h = function(d) {return(d)}
  } else {
    h = h_est
    h_d = h_d_est
  }

  if (length(C)==0) {
    data2 = cbind(Y,X,h(r))
  } else {
    data2 = cbind(Y,X,g_est(C),h(r))
  }
  names(data2)[ncol(data2)-1] = 'g'
  names(data2)[ncol(data2)] = 'h_r'
  # Stage II -- Estimation
  if (family == 'gaussian') {fa = gaussian()} else if (family == 'binary') {fa = binomial()}
  if (length(C) == 0) {
    c2  = mgcv::gam(y ~ s(x,bs="cr")+h_r, data = data2, family = fa)
  } else {
    c2  = mgcv::gam(y ~ s(x,bs="cr")+g+h_r, data = data2, family = fa)
  }
  mgcv::plot.gam(c2,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90",
           cex.lab = 1.3,cex.axis=1.3)

  if (family == 'gaussian') {
    if (action == 0) {
      V = var_cal_gaussian_r(c1,c2,data2$y)
    } else {
      V = var_cal_gaussian_h(c1,c2,data2$y,h_d)
    }
  } else if (family == 'binary') {
    if (action ==0 ) {
      V = var_cal_binomial_r(c1,c2,data2$y)
    } else {
      V = var_cal_binomial(c1,c2,data2$y,h_est = h)
    }
  }

  p = p_cal(c2,V,1)

  return(p)

}





















