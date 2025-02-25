rm(list=ls())
library(mgcv)
n1 = 5000; n2 = 10000; n3 = 20000; n4 = 50000



# 在此处更改样本量和工具变量强度
N = n4
PVE = 0.25
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
x = x+10
f = x     
y = 1+f+c+delta2

### nlmr
library(nlmr)
c = as.data.frame(c,ncol=1)
summary(x)

X = x
fp = fracpoly_mr(y, X, z, c, family = "gaussian", q = 10, d = "both", ci = "model_se", fig = TRUE)
summary(fp)
plm = piecewise_mr(y, X, z, c, family = "gaussian", q = 10, nboot = 50, fig = TRUE)
summary(plm)


set.seed(1016)
z = rnorm(N)
u = rnorm(N)
c = rnorm(N)
r1 = rnorm(N)
e = rnorm(N)
delta1 = u + r1
delta2 = delta1 + e
x = 1+bz*z+c+delta1
x = x+10
f = (x/3)^2  
y = 1+f+c+delta2

### nlmr
library(nlmr)
c = as.data.frame(c,ncol=1)
summary(x)

X = x
fp = fracpoly_mr(y, X, z, c, family = "gaussian", q = 10, d = "both", ci = "model_se", fig = TRUE)
summary(fp)
plm = piecewise_mr(y, X, z, c, family = "gaussian", q = 10, nboot = 50, fig = TRUE)
summary(plm)



set.seed(1016)
z = rnorm(N)
u = rnorm(N)
c = rnorm(N)
r1 = rnorm(N)
e = rnorm(N)
delta1 = u + r1
delta2 = delta1 + e
x = 1+bz*z+c+delta1
x = x+10
f = sin(x) 
y = 1+f+c+delta2

### nlmr
library(nlmr)
c = as.data.frame(c,ncol=1)
summary(x)

X = x
fp = fracpoly_mr(y, X, z, c, family = "gaussian", q = 10, d = "both", ci = "model_se", fig = TRUE)
summary(fp)
plm = piecewise_mr(y, X, z, c, family = "gaussian", q = 10, nboot = 50, fig = TRUE)
summary(plm)



set.seed(1016)
z = rnorm(N)
u = rnorm(N)
c = rnorm(N)
r1 = rnorm(N)
e = rnorm(N)
delta1 = u + r1
delta2 = delta1 + e
x = 1+bz*z+c+delta1
x = x+10
f = exp(x/3) 
y = 1+f+c+delta2

### nlmr
library(nlmr)
c = as.data.frame(c,ncol=1)
summary(x)

X = x
fp = fracpoly_mr(y, X, z, c, family = "gaussian", q = 10, d = "both", ci = "model_se", fig = TRUE)
summary(fp)
plm = piecewise_mr(y, X, z, c, family = "gaussian", q = 10, nboot = 50, fig = TRUE)
summary(plm)