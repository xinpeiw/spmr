
library(data.table)
a = fread('r_25_n1.txt')
apply(a,2,mean)
apply(a,2,sd)
a = fread('r_50_n1.txt')
apply(a,2,mean)
apply(a,2,sd)
a = fread('r_75_n1.txt')
apply(a,2,mean)
apply(a,2,sd)

a = fread('r_25_n2.txt')
apply(a,2,mean)
apply(a,2,sd)
a = fread('r_50_n2.txt')
apply(a,2,mean)
apply(a,2,sd)
a = fread('r_75_n2.txt')
apply(a,2,mean)
apply(a,2,sd)


a = fread('r_25_n3.txt')
apply(a,2,mean)
apply(a,2,sd)
a = fread('r_50_n3.txt')
apply(a,2,mean)
apply(a,2,sd)
a = fread('r_75_n3.txt')
apply(a,2,mean)
apply(a,2,sd)


a = fread('r_25_n4.txt')
apply(a,2,mean)
apply(a,2,sd)
a = fread('r_50_n4.txt')
apply(a,2,mean)
apply(a,2,sd)
a = fread('r_75_n4.txt')
apply(a,2,mean)
apply(a,2,sd)





## 2sp 
n1 = 5000; n2 = 10000; n3 = 20000; n4 = 50000

var_est_sp = function(N,PVE) {
  bz = sqrt(3*PVE/(1-PVE))
  set.seed(1016)
  z = rnorm(N)
  u = rnorm(N)
  c = rnorm(N)
  r1 = rnorm(N)
  e = rnorm(N)
  
  delta1 = u + r1
  delta2 = delta1 + e
  
  x = 1+bz*z+c+delta1
  f = exp(x/3)  
  y = 1+f+c+delta2
  
  s1 = lm(f ~ z+c)
  p = s1$"fitted.values"
  
  s2 = lm(y ~ p+c)
  
  V = matrix(c(rep(1,N),p,c),ncol=3)
  return(sqrt((solve(t(V)%*%V) * var(delta2))[2,2]))
}
var_est_sp(5000,0.25)
var_est_sp(5000,0.50)
var_est_sp(5000,0.75)

var_est_sp(10000,0.25)
var_est_sp(10000,0.50)
var_est_sp(10000,0.75)

var_est_sp(20000,0.25)
var_est_sp(20000,0.50)
var_est_sp(20000,0.75)

var_est_sp(50000,0.25)
var_est_sp(50000,0.50)
var_est_sp(50000,0.75)



## CF
var_est_cf = function(N,PVE){
  bz = sqrt(3*PVE/(1-PVE))
  set.seed(1016)
  z = rnorm(N)
  u = rnorm(N)
  c = rnorm(N)
  r1 = rnorm(N)
  e = rnorm(N)
  
  delta1 = u + r1
  delta2 = delta1 + e
  
  x = 1+bz*z+c+delta1
  f = exp(x/3)
  y = 1+f+c+delta2
  
  c1 = lm(x ~ z+c)
  r = c1$residuals
  
  c2 = lm(y ~ f+c+r)
  
  W = matrix(cbind(rep(1,N),f,c,r),ncol = 4)
  rr = matrix(delta1-r,ncol=1)
  ehat = c2$residuals
  r2 = matrix(ehat,ncol=1)
  dbeta0 = as.numeric(var(delta1)/N)
  dbetaz = as.numeric(var(delta1)/N/var(z))
  dbetac = as.numeric(var(delta1)/N/var(c))
  Z = matrix(z,ncol=1)
  C = matrix(c,ncol=1)
  I = matrix(rep(1,N),ncol=1)
  rho = 1
  
  v = (solve(t(W)%*%W) * (as.numeric(var(r2)) + 
                             rho^2 * solve(t(W)%*%W) %*% (t(W)%*%I%*%t(I)%*%W) * dbeta0 +
                             rho^2 * solve(t(W)%*%W) %*% (t(W)%*%Z%*%t(Z)%*%W) * dbetaz + 
                             rho^2 * solve(t(W)%*%W) %*% (t(W)%*%C%*%t(C)%*%W) * dbetac)[2,2])[2,2]
  return(sqrt(v))
}

var_est_cf(5000,0.25)
var_est_cf(5000,0.50)
var_est_cf(5000,0.75)

var_est_cf(10000,0.25)
var_est_cf(10000,0.50)
var_est_cf(10000,0.75)

var_est_cf(20000,0.25)
var_est_cf(20000,0.50)
var_est_cf(20000,0.75)

var_est_cf(50000,0.25)
var_est_cf(50000,0.50)
var_est_cf(50000,0.75)

