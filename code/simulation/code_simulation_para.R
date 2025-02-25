n1 = 5000; n2 = 10000; n3 = 20000; n4 = 50000

simu = function(PVE,N) {
  res_CF = c()
  res_SP = c()
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
    f = exp(x/3)       # 在此处更改f的形式
    y = 1+f+c+delta2
    
    # CF
    r = lm(x ~ z+c)$residuals
    res_CF[i] = lm(y ~ f+c+r)$coefficients[[2]]
    
    p = lm(f ~ z+c)$"fitted.values"
    res_SP[i] = lm(y ~ p+c)$coefficients[[2]]
  }
  return(cbind(res_CF,res_SP))
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

setwd('C:/ShareCache/王鑫培_2111110203/AAA小论文/2 - nonlinearMR/Simulation/exp')
# 根据f的形式在此处更改文件夹名称
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