n1 = 5000; n2 = 10000; n3 = 20000; n4 = 50000

simu = function(PVE,N) {
  res_CF = c()
  res_false = c()
  res_non = c()
  bz = sqrt(3*PVE/(1-PVE))
  set.seed(1016)
  for (i in 1:1000) {
    if (i %% 20 == 0) { print(i) }
    z = rnorm(N)
    u = rnorm(N)
    c = rnorm(N)
    r1 = rnorm(N)
    e = rnorm(N)
    
    delta1 = u + r1
    trans = function(data){
      return(exp(data/3))
    }
    h = trans(delta1)
    delta2 = h + e
    
    x = 1+bz*z+c+delta1
    f = sin(x)
    y = 1+f+c+delta2
    
    # CF
    r = lm(x ~ z+c)$residuals
    res_CF[i] = lm(y ~ f+c+I(trans(r)))$coefficients[[2]]
    
    # false
    res_false[i] = lm(y ~ f+c+r)$coefficients[[2]]
    
    # non-para
    c1 = lm(x ~ z+c)
    r = c1$residuals
    V = model.matrix(c1)
    omega = vcov(c1)
    M = 500
    hat_h = c()
    vr = rowSums((V%*%omega)*V)
    
    d = cbind(r,vr)
    samp = function(data){
      sam = rnorm(M,data[[1]],sqrt(data[[2]]))
      return(mean(trans(sam)))
    }
    hat_h = apply(d, 1, samp)
    
    res_non[i] = lm(y ~ f+c+hat_h)$coefficients[[2]]
    
  }
  return(cbind(res_CF,res_false,res_non))
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


setwd('C:/ShareCache/王鑫培_2111110203/AAA小论文/4 - CF_nolinear/simulation/exp')
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

r_25_n1 = fread('./exp/r_25_n1.txt')
r_25_n1 = r_25_n1[,c(1,3,2)]
names(r_25_n1) = c('参数估计','非参数估计','线性估计')
boxplot(r_25_n1,
        xlab = "估计方法",ylab = "系数估计值" )


## spline
setwd("C:/ShareCache/王鑫培_2111110203/AAA小论文/4 - CF_nolinear/simulation/附录表")
n1 = 5000; n2 = 10000; n3 = 20000; n4 = 50000
PVE = 0.25
N = n1
bz = sqrt(3*PVE/(1-PVE))
set.seed(1016)
z = rnorm(N)
u = rnorm(N)
c = rnorm(N)
r1 = rnorm(N)
e = rnorm(N)

delta1 = u + r1
trans = function(data){
  return(exp(data/3))
}
h = trans(delta1)
delta2 = h + e

x = 1+bz*z+c+delta1
f = sin(x)
y = 1+f+c+delta2


residuals = lm(x ~ z+c)$residuals
library(mgcv)
X = x
o = gam(y ~ s(X,bs="cr")+c+s(residuals,bs="cr"),family = "gaussian")
plot.gam(o,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90",
         cex.lab = 1.3,cex.axis=1.3)











