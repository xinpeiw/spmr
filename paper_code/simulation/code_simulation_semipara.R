rm(list=ls())
library(mgcv)
n1 = 5000; n2 = 10000; n3 = 20000; n4 = 50000


N = n1
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
x = 10+bz*z+c+delta1
f = x
y = 1+f+c+delta2

r = lm(x ~ z+c)$residuals
X = x
o_x = gam(y ~ s(X,bs="cr")+c+r, family = gaussian())
plot.gam(o_x,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         cex.lab = 1.3,cex.axis=1.3)


N = n1
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
x = 10+bz*z+c+delta1
f = (x/3)^2
y = 1+f+c+delta2

r = lm(x ~ z+c)$residuals
X = x
o_x2 = gam(y ~ s(X,bs="cr")+c+r, family = gaussian())
plot.gam(o_x2,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         cex.lab = 1.3,cex.axis=1.3)


N = n1
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
x = 10+bz*z+c+delta1
f = sin(x)
y = 1+f+c+delta2

r = lm(x ~ z+c)$residuals
X = x
o_sin = gam(y ~ s(X,bs="cr")+c+r, family = gaussian())
plot.gam(o_sin,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         cex.lab = 1.3,cex.axis=1.3)



N = n1
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
x = 10+bz*z+c+delta1
f = exp(x/3)
y = 1+f+c+delta2

r = lm(x ~ z+c)$residuals
X = x
o_exp = gam(y ~ s(X,bs="cr")+c+r, family = gaussian())
plot.gam(o_exp,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         cex.lab = 1.3,cex.axis=1.3)





par(mfrow=c(2,2))
plot.gam(o_x,rug = TRUE,se=TRUE,n=N,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         family = "serif")
title(main = expression('f'*'(X)'*' '*'='*' '*X), family = "serif", cex.main = 1, font.main= 5)

plot.gam(o_x2,rug = TRUE,se=TRUE,n=N,jit = TRUE,shade = TRUE,shade.col = "grey90", family = "serif")
title(main = expression('f'*'(X)'*' '*'='*' '*(X/3)^2), family = "serif", cex.main = 1, font.main= 5)

plot.gam(o_sin,rug = TRUE,se=TRUE,n=N,jit = TRUE,shade = TRUE,shade.col = "grey90", family = "serif")
title(main = expression('f'*'(X)'*' '*'='*' '*'sin'*'(X)'), family = "serif", cex.main = 1, font.main= 5)

plot.gam(o_exp,rug = TRUE,se=TRUE,n=N,jit = TRUE,shade = TRUE,shade.col = "grey90", family = "serif")
title(main = expression('f'*'(X)'*' '*'='*' '*e^{X/3}), family = "serif", cex.main = 1, font.main= 5)

