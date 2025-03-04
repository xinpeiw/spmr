### Invalid IV ###

# uncorrelated horizontal pleiotropy #
rm(list=ls())
n1 = 5000; n2 = 10000; n3 = 20000; n4 = 50000

simu = function(PVE,N) {
  res_CF = c()
  bz = sqrt(3*PVE/(1-PVE))
  set.seed(1016)
  for (i in 1:1000) {
    z = rnorm(N)
    u = rnorm(N)
    c = rnorm(N)
    r1 = rnorm(N)
    e = rnorm(N)
    
    delta1 = u + r1
    delta2 = delta1 + e
    
    x = 1+bz*z+c+delta1
    f = exp(x/3)
    y = 1+f+c+delta2+z
    
    # CF
    r = lm(x ~ z+c)$residuals
    res_CF[i] = lm(y ~ f+c+r+z)$coefficients[[2]]
  }
  return(res_CF)
}
r_25_n1 = as.data.frame(simu(0.25,n1))
r_25_n2 = as.data.frame(simu(0.25,n2))
r_25_n3 = as.data.frame(simu(0.25,n3))
r_25_n4 = as.data.frame(simu(0.25,n4))

r_50_n1 = as.data.frame(simu(0.50,n1))
r_50_n2 = as.data.frame(simu(0.50,n2))
r_50_n3 = as.data.frame(simu(0.50,n3))
r_50_n4 = as.data.frame(simu(0.50,n4))

r_75_n1 = as.data.frame(simu(0.75,n1))
r_75_n2 = as.data.frame(simu(0.75,n2))
r_75_n3 = as.data.frame(simu(0.75,n3))
r_75_n4 = as.data.frame(simu(0.75,n4))

write.table(r_25_n1,'r_25_n1.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_25_n2,'r_25_n2.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_25_n3,'r_25_n3.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_25_n4,'r_25_n4.txt',row.names = FALSE, sep=' ', quote = FALSE)

write.table(r_50_n1,'r_50_n1.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_50_n2,'r_50_n2.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_50_n3,'r_50_n3.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_50_n4,'r_50_n4.txt',row.names = FALSE, sep=' ', quote = FALSE)

write.table(r_75_n1,'r_75_n1.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_75_n2,'r_75_n2.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_75_n3,'r_75_n3.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_75_n4,'r_75_n4.txt',row.names = FALSE, sep=' ', quote = FALSE)

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




## Variance estimation ##

var_est = function(PVE,N) {
  set.seed(1016)
  bz = sqrt(3*PVE/(1-PVE))
  
  z = rnorm(N)
  u = rnorm(N)
  c = rnorm(N)
  r1 = rnorm(N)
  e = rnorm(N)
  
  delta1 = u + r1
  delta2 = delta1 + e
  
  x = 1+bz*z+c+delta1
  f = exp(x/3)
  y = 1+f+c+delta2+z
  
  c1 = lm(x ~ z+c)
  r = c1$residuals
  c2 = lm(y ~ f+c+r+z)
  
  W = as.matrix(cbind(rep(1,N),f,c,r,z),ncol=5)
  I = matrix(rep(1,N),ncol=1)
  Z = matrix(z,ncol=1)
  C = matrix(c,ncol=1)
  rho = 1
  
  dbeta0 = vcov(c1)[1,1]
  dbetaz = vcov(c1)[2,2]
  dbetac = vcov(c1)[3,3]
  dbeta0z = vcov(c1)[1,2]
  dbeta0c = vcov(c1)[1,3]
  dbetazc = vcov(c1)[2,3]
  
  
  ehat = c2$residuals
  
  v = (solve(t(W)%*%W) * (as.numeric(var(ehat)) + 
                        rho^2 * solve(t(W)%*%W) %*% (t(W)%*%I%*%t(I)%*%W) * dbeta0 +
                        rho^2 * solve(t(W)%*%W) %*% (t(W)%*%Z%*%t(Z)%*%W) * dbetaz + 
                        rho^2 * solve(t(W)%*%W) %*% (t(W)%*%C%*%t(C)%*%W) * dbetac + 
                        2*rho^2 * solve(t(W)%*%W) %*% (t(W)%*%I%*%t(C)%*%W) * dbeta0c + 
                        2*rho^2 * solve(t(W)%*%W) %*% (t(W)%*%I%*%t(Z)%*%W) * dbeta0z + 
                        2*rho^2 * solve(t(W)%*%W) %*% (t(W)%*%Z%*%t(C)%*%W) * dbetazc
  )[2,2])[2,2]
  
  return(sqrt(v))
  
}

var_est(0.25,5000)
var_est(0.50,5000)
var_est(0.75,5000)

var_est(0.25,10000)
var_est(0.50,10000)
var_est(0.75,10000)

var_est(0.25,20000)
var_est(0.50,20000)
var_est(0.75,20000)

var_est(0.25,50000)
var_est(0.50,50000)
var_est(0.75,50000)








# Correlated horizontal pleiotropy #
rm(list=ls())
n1 = 5000; n2 = 10000; n3 = 20000; n4 = 50000

simu = function(PVE,N) {
  res_CF = c()
  bz = sqrt(3*PVE/(1-PVE))
  set.seed(1016)
  for (i in 1:1000) {
    z = rnorm(N)
    u = z + rnorm(N)
    c = rnorm(N)
    r1 = rnorm(N)
    e = rnorm(N)
    
    delta1 = u + r1
    delta2 = delta1 + e
    
    x = 1+bz*z+c+delta1
    f = exp(x/3)
    y = 1+f+c+delta2
    
    # CF
    r = lm(x ~ z+c)$residuals
    res_CF[i] = lm(y ~ f+c+r+z)$coefficients[[2]]
  }
  return(res_CF)
}
r_25_n1 = as.data.frame(simu(0.25,n1))
r_25_n2 = as.data.frame(simu(0.25,n2))
r_25_n3 = as.data.frame(simu(0.25,n3))
r_25_n4 = as.data.frame(simu(0.25,n4))

r_50_n1 = as.data.frame(simu(0.50,n1))
r_50_n2 = as.data.frame(simu(0.50,n2))
r_50_n3 = as.data.frame(simu(0.50,n3))
r_50_n4 = as.data.frame(simu(0.50,n4))

r_75_n1 = as.data.frame(simu(0.75,n1))
r_75_n2 = as.data.frame(simu(0.75,n2))
r_75_n3 = as.data.frame(simu(0.75,n3))
r_75_n4 = as.data.frame(simu(0.75,n4))


write.table(r_25_n1,'r_25_n1.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_25_n2,'r_25_n2.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_25_n3,'r_25_n3.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_25_n4,'r_25_n4.txt',row.names = FALSE, sep=' ', quote = FALSE)

write.table(r_50_n1,'r_50_n1.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_50_n2,'r_50_n2.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_50_n3,'r_50_n3.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_50_n4,'r_50_n4.txt',row.names = FALSE, sep=' ', quote = FALSE)

write.table(r_75_n1,'r_75_n1.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_75_n2,'r_75_n2.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_75_n3,'r_75_n3.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_75_n4,'r_75_n4.txt',row.names = FALSE, sep=' ', quote = FALSE)


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


## Variance estimation ##

var_est = function(PVE,N) {
  set.seed(1016)
  bz = sqrt(3*PVE/(1-PVE))
  
  z = rnorm(N)
  u = z + rnorm(N)
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
  c2 = lm(y ~ f+c+r+z)
  
  W = as.matrix(cbind(rep(1,N),f,c,r,z),ncol=5)
  I = matrix(rep(1,N),ncol=1)
  Z = matrix(z,ncol=1)
  C = matrix(c,ncol=1)
  rho = 1
  
  dbeta0 = vcov(c1)[1,1]
  dbetaz = vcov(c1)[2,2]
  dbetac = vcov(c1)[3,3]
  dbeta0z = vcov(c1)[1,2]
  dbeta0c = vcov(c1)[1,3]
  dbetazc = vcov(c1)[2,3]
  
  
  ehat = c2$residuals
  
  v = (solve(t(W)%*%W) * (as.numeric(var(ehat)) + 
                            rho^2 * solve(t(W)%*%W) %*% (t(W)%*%I%*%t(I)%*%W) * dbeta0 +
                            rho^2 * solve(t(W)%*%W) %*% (t(W)%*%Z%*%t(Z)%*%W) * dbetaz + 
                            rho^2 * solve(t(W)%*%W) %*% (t(W)%*%C%*%t(C)%*%W) * dbetac + 
                            2*rho^2 * solve(t(W)%*%W) %*% (t(W)%*%I%*%t(C)%*%W) * dbeta0c + 
                            2*rho^2 * solve(t(W)%*%W) %*% (t(W)%*%I%*%t(Z)%*%W) * dbeta0z + 
                            2*rho^2 * solve(t(W)%*%W) %*% (t(W)%*%Z%*%t(C)%*%W) * dbetazc
  )[2,2])[2,2]
  
  return(sqrt(v))
  
}

var_est(0.25,5000)
var_est(0.50,5000)
var_est(0.75,5000)

var_est(0.25,10000)
var_est(0.50,10000)
var_est(0.75,10000)

var_est(0.25,20000)
var_est(0.50,20000)
var_est(0.75,20000)

var_est(0.25,50000)
var_est(0.50,50000)
var_est(0.75,50000)












# Both Uncorrelated and Correlated horizontal pleiotropy #
rm(list=ls())
n1 = 5000; n2 = 10000; n3 = 20000; n4 = 50000

simu = function(PVE,N) {
  res_CF = c()
  bz = sqrt(3*PVE/(1-PVE))
  set.seed(1016)
  for (i in 1:1000) {
    z = rnorm(N)
    u = z + rnorm(N)
    c = rnorm(N)
    r1 = rnorm(N)
    e = rnorm(N)
    
    delta1 = u + r1
    delta2 = delta1 + e
    
    x = 1+bz*z+c+delta1
    f = exp(x/3)
    y = 1+f+c+delta2+z
    
    # CF
    r = lm(x ~ z+c)$residuals
    res_CF[i] = lm(y ~ f+c+r+z)$coefficients[[2]]
  }
  return(res_CF)
}
r_25_n1 = as.data.frame(simu(0.25,n1))
r_25_n2 = as.data.frame(simu(0.25,n2))
r_25_n3 = as.data.frame(simu(0.25,n3))
r_25_n4 = as.data.frame(simu(0.25,n4))

r_50_n1 = as.data.frame(simu(0.50,n1))
r_50_n2 = as.data.frame(simu(0.50,n2))
r_50_n3 = as.data.frame(simu(0.50,n3))
r_50_n4 = as.data.frame(simu(0.50,n4))

r_75_n1 = as.data.frame(simu(0.75,n1))
r_75_n2 = as.data.frame(simu(0.75,n2))
r_75_n3 = as.data.frame(simu(0.75,n3))
r_75_n4 = as.data.frame(simu(0.75,n4))


write.table(r_25_n1,'r_25_n1.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_25_n2,'r_25_n2.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_25_n3,'r_25_n3.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_25_n4,'r_25_n4.txt',row.names = FALSE, sep=' ', quote = FALSE)

write.table(r_50_n1,'r_50_n1.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_50_n2,'r_50_n2.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_50_n3,'r_50_n3.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_50_n4,'r_50_n4.txt',row.names = FALSE, sep=' ', quote = FALSE)

write.table(r_75_n1,'r_75_n1.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_75_n2,'r_75_n2.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_75_n3,'r_75_n3.txt',row.names = FALSE, sep=' ', quote = FALSE)
write.table(r_75_n4,'r_75_n4.txt',row.names = FALSE, sep=' ', quote = FALSE)



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


## Variance estimation ##

var_est = function(PVE,N) {
  set.seed(1016)
  bz = sqrt(3*PVE/(1-PVE))
  
  z = rnorm(N)
  u = z + rnorm(N)
  c = rnorm(N)
  r1 = rnorm(N)
  e = rnorm(N)
  
  delta1 = u + r1
  delta2 = delta1 + e
  
  x = 1+bz*z+c+delta1
  f = exp(x/3)
  y = 1+f+c+delta2+z
  
  c1 = lm(x ~ z+c)
  r = c1$residuals
  c2 = lm(y ~ f+c+r+z)
  
  W = as.matrix(cbind(rep(1,N),f,c,r,z),ncol=5)
  I = matrix(rep(1,N),ncol=1)
  Z = matrix(z,ncol=1)
  C = matrix(c,ncol=1)
  rho = 1
  
  dbeta0 = vcov(c1)[1,1]
  dbetaz = vcov(c1)[2,2]
  dbetac = vcov(c1)[3,3]
  dbeta0z = vcov(c1)[1,2]
  dbeta0c = vcov(c1)[1,3]
  dbetazc = vcov(c1)[2,3]
  
  
  ehat = c2$residuals
  
  v = (solve(t(W)%*%W) * (as.numeric(var(ehat)) + 
                            rho^2 * solve(t(W)%*%W) %*% (t(W)%*%I%*%t(I)%*%W) * dbeta0 +
                            rho^2 * solve(t(W)%*%W) %*% (t(W)%*%Z%*%t(Z)%*%W) * dbetaz + 
                            rho^2 * solve(t(W)%*%W) %*% (t(W)%*%C%*%t(C)%*%W) * dbetac + 
                            2*rho^2 * solve(t(W)%*%W) %*% (t(W)%*%I%*%t(C)%*%W) * dbeta0c + 
                            2*rho^2 * solve(t(W)%*%W) %*% (t(W)%*%I%*%t(Z)%*%W) * dbeta0z + 
                            2*rho^2 * solve(t(W)%*%W) %*% (t(W)%*%Z%*%t(C)%*%W) * dbetazc
  )[2,2])[2,2]
  
  return(sqrt(v))
  
}

var_est(0.25,5000)
var_est(0.50,5000)
var_est(0.75,5000)

var_est(0.25,10000)
var_est(0.50,10000)
var_est(0.75,10000)

var_est(0.25,20000)
var_est(0.50,20000)
var_est(0.75,20000)

var_est(0.25,50000)
var_est(0.50,50000)
var_est(0.75,50000)








