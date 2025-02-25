
# 所需函数
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


var_cal_binomial = function(object1,object2,y,i) 
  # object1: stage 1 model
  # object2: gam model
  # i: Number of the interested smoothing term
{
  Xp = model.matrix(object2)
  c2 = glm(y ~ Xp-1, family = "binomial",maxit=50)
  B = matrix(coef(c2),ncol=1)
  etahat = c2$linear.predictors
  miuhat = c2$fitted.values
  gz = Xp * as.numeric(y-miuhat)
  
  V = model.matrix(object1)
  r = object1$residuals
  N = length(y)
  fiez = as.matrix(V * r, ncol=ncol(V)) %*% solve(t(V)%*%V) * N
  a = y - miuhat - etahat*miuhat*(1-miuhat)
  index.r = which(colnames(Xp)=='r')
  Ggamma = matrix(c(rep(0,ncol(V)*(index.r-1)),-apply(as.numeric(a)*V,2,mean),rep(0,ncol(V)*(ncol(Xp)-index.r)))
                  ,ncol=ncol(Xp))
  A = gz + fiez %*% Ggamma
  
  S = (object2$smooth[[1]])$S[[1]]
  S = rbind(matrix(rep(0,ncol(Xp)*(ncol(Xp)-nrow(S))),nrow = (ncol(Xp)-nrow(S))),
            cbind(matrix(rep(0,(ncol(Xp)-nrow(S))*nrow(S)),ncol = (ncol(Xp)-nrow(S))),S))
  lambda = object2$sp[[1]]
  W = Xp
  V = vcov(object2) %*% crossprod(A) %*% solve(solve(vcov(object2))-lambda*S) 
  
  return(V)
}




# P_value calculation
p_cal = function(object,V,i) {  # object: gam; V: variance; i: No.smoothing term
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



## Binomial simulation
set.seed(1016)
p1 = c()
p2 = c()
N = 5000
library(mgcv)
n1 = 5000; n2 = 10000; n3 = 20000; n4 = 50000
p = c()
p_l = c()
simu = function(PVE,N) {
  bz = sqrt(3*PVE/(1-PVE))
  set.seed(1016)
  for (i in 1:500) {
    # data generating
    z = rnorm(N)
    u = rnorm(N)
    c = rnorm(N)
    r1 = rnorm(N)
    e = rnorm(N)
    
    delta1 = u + r1
    delta2 = delta1 + e
    
    x = 1+bz*z+c+delta1
    f = 0.5*(x/3)^2
    eta = 1+f+c+delta1
    miu = exp(eta)/(1+exp(eta))
    y = rbinom(N,1,miu)
    
    # model fitting
    c1 = lm(x ~ z + c)
    r = c1$residuals
    o = gam(y ~ s(x,bs="cr") + c +r, family = "binomial")
    
    V1 = var_cal_binomial(c1,o,y,1)
    p[i] = p_cal(o,V1,1)
    
    
    # 线性假设下
    c2 = glm(y ~ x+c+r,family = "binomial",maxit=50)
    W = model.matrix(c2)
    rho = coef(c2)[[4]]
    r2 = matrix(c2$residuals,ncol=1)
    V = model.matrix(c1)
    
    etahat = c2$linear.predictors
    miuhat = c2$fitted.values
    gz = W * as.numeric(y-miuhat)
    fiez = as.matrix(V * r, ncol=ncol(V)) %*% solve(crossprod(V)) * N
    a = y - miuhat - etahat*miuhat*(1-miuhat)
    index.r = which(colnames(W)=='r')
    Ggamma = matrix(c(rep(0,ncol(V)*(index.r-1)),-apply(as.numeric(a)*V,2,mean),rep(0,ncol(V)*(ncol(W)-index.r)))
      ,ncol=ncol(W))
    A = gz + fiez %*% Ggamma
    Wei = diag(c2$weights)
    se = sqrt((vcov(c2) %*% crossprod(A) %*% vcov(c2))[2,2])
    p_l[i] = pnorm(-abs(coef(c2)[[2]]/se))*2
    print(i)
  }
  return(cbind(p,p_l))
}

r_25_n1_x2 = simu(0.25,n1)
r_25_n1_x2 = as.data.frame(r_25_n1_x2)
summary(r_25_n1_x2[,1])
summary(r_25_n1_x2[,2])
p_u = runif(1000)
qqplot(p_u,r_25_n1_x2[,1],type='l',lty=1,ylim = c(0,1),
       xlab = "Uniform quantiles",ylab="Ordered p-values",cex.lab = 1.5,cex.axis=1.5)
par(new=T)
qqplot(p_u,r_25_n1_x2[,2],type='l',lty=2,ylim = c(0,1),
       xlab = "Uniform quantiles",ylab="Ordered p-values",cex.lab = 1.5,cex.axis=1.5)


## f = 0
library(mgcv)
n1 = 5000; n2 = 10000; n3 = 20000; n4 = 50000
simu = function(PVE,N) {
  bz = sqrt(3*PVE/(1-PVE))
  set.seed(1016)
  p = c()
  for (i in 1:500) {
    # data generating
    print(i)
    z = rnorm(N)
    u = rnorm(N)
    c = rnorm(N)
    r1 = rnorm(N)
    e = rnorm(N)
    
    delta1 = u + r1
    delta2 = delta1 + e
    
    x = 1+bz*z+c+delta1
    f = 0
    eta = 1+f+c+delta1
    miu = exp(eta)/(1+exp(eta))
    y = rbinom(N,1,miu)
    
    # model fitting
    c1 = lm(x ~ z + c)
    r = c1$residuals
    o = gam(y ~ s(x,bs="cr") + c +r, family = "binomial")
    
    V1 = var_cal_binomial(c1,o,y,1)
    p[i] = p_cal(o,V1,1)
  }
  return(p)
}


r_25_n1 = simu(0.25,n1)
r_25_n2 = simu(0.25,n2)
r_25_n3 = simu(0.25,n3)
r_25_n4 = simu(0.25,n4)

r_50_n1 = simu(0.50,n1)
r_50_n2 = simu(0.50,n2)
r_50_n3 = simu(0.50,n3)
r_50_n4 = simu(0.50,n4)

r_75_n1 = simu(0.75,n1)
r_75_n2 = simu(0.75,n2)
r_75_n3 = simu(0.75,n3)
r_75_n4 = simu(0.75,n4)


r_25 = cbind(r_25_n1,r_25_n2,r_25_n3,r_25_n4)
write.table(r_25,'C:/ShareCache/王鑫培_2111110203/AAA小论文/5 - logistic/simulation/smoothingP0_0.25.txt',
            sep=' ',quote = FALSE,row.names = FALSE)

r_75 = cbind(r_75_n1,r_75_n2,r_75_n3,r_75_n4)
write.table(r_75,
            'C:/ShareCache/王鑫培_2111110203/AAA小论文/5 - logistic/simulation/smoothingP0_0.75.txt',
            sep=' ',quote = FALSE,row.names = FALSE)

r_50 = cbind(r_50_n1,r_50_n2,r_50_n3,r_50_n4)
write.table(r_50,
            'C:/ShareCache/王鑫培_2111110203/AAA小论文/5 - logistic/simulation/smoothingP0_0.50.txt',
            sep=' ',quote = FALSE,row.names = FALSE)

library(data.table)
r_25 = fread('C:/ShareCache/王鑫培_2111110203/AAA小论文/5 - logistic/simulation/smoothingP0_0.25.txt')
r_50 = fread('C:/ShareCache/王鑫培_2111110203/AAA小论文/5 - logistic/simulation/smoothingP0_0.50.txt')
r_75 = fread('C:/ShareCache/王鑫培_2111110203/AAA小论文/5 - logistic/simulation/smoothingP0_0.75.txt')
r = cbind(r_25,r_50,r_75)
alpha = function(data){
  p = length(data[which(data<0.05)]) / length(data)
}


# 画柱状图
count = as.numeric(apply(r, 2, alpha))
names = c('25%, 5000','25%,10000','25%,20000','25%,50000',
          '50%, 5000','50%,10000','50%,20000','50%,50000',
          '75%, 5000','75%,10000','75%,20000','75%,50000')
par(omi=c(0.1,0,0,0),mar=c(9.5,5,4,3))
k = barplot(p,
            xlab = " ",ylab = "I类错误率",
            names.arg = names, xaxt="n", space=1)
text(x=k,y=-0.003,srt = 45, adj= 1, xpd = TRUE, labels = names)
text(x=13.5, y=-0.02,xpd=TRUE, labels="模拟设置")