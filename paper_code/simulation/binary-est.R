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
    
    x = 1+bz*z+c+delta1
    f = exp(x/3)
    eta = 1+f+c+delta1
    miu = exp(eta)/(1+exp(eta))
    y = rbinom(N,1,miu)
    
    # CF
    r = lm(x ~ z+c)$residuals
    res_CF[i] = glm(y ~ f+c+r, family = "binomial",maxit=50)$coefficients[[2]]
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

mean(r_25_n1[,1]);sd(r_25_n1[,1]);mean(r_50_n1[,1]);sd(r_50_n1[,1]);mean(r_75_n1[,1]);sd(r_75_n1[,1])
mean(r_25_n2[,1]);sd(r_25_n2[,1]);mean(r_50_n2[,1]);sd(r_50_n2[,1]);mean(r_75_n2[,1]);sd(r_75_n2[,1])
mean(r_25_n3[,1]);sd(r_25_n3[,1]);mean(r_50_n3[,1]);sd(r_50_n3[,1]);mean(r_75_n3[,1]);sd(r_75_n3[,1])
mean(r_25_n4[,1]);sd(r_25_n4[,1]);mean(r_50_n4[,1]);sd(r_50_n4[,1]);mean(r_75_n4[,1]);sd(r_75_n4[,1])



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






var_est = function(PVE,N) {
  bz = sqrt(3*PVE/(1-PVE))
  set.seed(1016)
  z = rnorm(N)
  u = rnorm(N)
  c = rnorm(N)
  r1 = rnorm(N)
  e = rnorm(N)
  
  delta1 = u + r1
  
  x = 1+bz*z+c+delta1
  f = exp(x/3)
  eta = 1+f+c+delta1
  miu = exp(eta)/(1+exp(eta))
  y = rbinom(N,1,miu)
  
  # CF
  c1 = lm(x ~ z+c)
  r = c1$residuals
  c2 = glm(y ~ f+c+r, family = "binomial",maxit=50)
  
  
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
  #Wei = diag(c2$weights)
  v = vcov(c2)
  se = sqrt((v %*% crossprod(A) %*% v)[2,2])
  
  
  
  return(se)
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



n1 = 5000; n2 = 10000; n3 = 20000; n4 = 50000
PVE = 0.75
N = n4
bz = sqrt(3*PVE/(1-PVE))
set.seed(1016)
z = rnorm(N)
u = rnorm(N)
c = rnorm(N)
r1 = rnorm(N)
e = rnorm(N)

delta1 = u + r1

x = 1+bz*z+c+delta1
f = x
eta = 1+f+c+delta1
miu = exp(eta)/(1+exp(eta))
y = rbinom(N,1,miu)

# CF
c1 = lm(x ~ z+c)
r = c1$residuals
X = x
c2 = gam(y ~ s(X,bs="cr")+c+r, family = binomial,maxit=50)
plot.gam(c2,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         cex.lab = 1.3,cex.axis=1.3)

