## 匹配病例-对照
explist = c('BMI','WC','HIP','WHR','BF','BMR','DBP','SBP','sleep_duration','alcohol','coffee')
outcomes = c('I11','I20','I21','I70','I251','I48','I50','I61','I63')
library(MatchIt)
library(data.table)
for (expi in 1:length(explist)) {
  exp = explist[expi]
  print(exp)
  exppath = paste0(paste0('data_',exp),'.txt')
  data = fread(exppath)
  data = as.data.frame(data)
  
  data$sex = as.factor(data$sex)
  data$centre = as.factor(data$centre)
  data$batch = as.factor(data$batch)
  data[,5:13]<-lapply(data[,5:13],as.factor)
  
  for (outi in 1:length(outcomes)) {
    data$out = data[,which(colnames(data)==outcomes[outi])]
    print(outcomes[outi])
    matchlist <- matchit(out ~ sex+age+batch+centre, data=data,
                         method   = "nearest",
                         distance = "glm",      #
                         caliper  = 0.05,       # 卡钳值
                         ratio    = 2,          # 1:N 匹配
                         replace  = F)          # 不替换
    
    matchdata<- match.data(matchlist,
                           group = "all",
                           distance = "distance",
                           weights = "weights",
                           subclass = "subclass",
                           data = NULL,
                           include.s.weights = TRUE,
                           drop.unmatched = TRUE)
    print(table(duplicated(matchdata$n_eid)))
    print(table(matchdata$out)) # 匹配后
    print(table(data$out))      # 匹配前
    
    a=table(matchdata$subclass);print(table(a))
    
    print(dim(matchdata))
    outpath = paste0(paste0(paste0('./matched/matched.data/',exp),paste0('.',outcomes[outi])),'.txt')
    write.table(matchdata,outpath,row.names = FALSE,sep=' ',quote = FALSE)
    
  }
}
write.table(matchdata,'./matched/matched.data/BMR.I70.1vs1.txt',row.names = FALSE,sep=' ',quote = FALSE)
## coffee
expi = 11
exp = explist[expi]
print(exp)
exppath = paste0(paste0('data_',exp),'.txt')
data = fread(exppath)
data = as.data.frame(data)

data$sex = as.factor(data$sex)
data$centre = as.factor(data$centre)
data$batch = as.factor(data$batch)
data[,5:14]<-lapply(data[,5:14],as.factor)

d = data[-which(data$coffee == -1 | data$coffee == -3),]
d1 = d[which(data$coffee == -10),]
d1 = d1[!is.na(d1$n_eid),]
d1$coffee = 0.5
d2 = d[-which(data$coffee == -10),]
data = rbind(d1,d2)
summary(data$coffee)

library(MatchIt)
for (outi in 1:length(outcomes)) {
  data$out = data[,which(colnames(data)==outcomes[outi])]
  print(outcomes[outi])
  matchlist <- matchit(out ~ sex+age+batch+centre, data=data,
                       method   = "nearest",
                       distance = "glm",      #
                       caliper  = 0.05,       # 卡钳值
                       ratio    = 2,          # 1:N 匹配
                       replace  = F)          # 不替换
  
  matchdata<- match.data(matchlist,
                         group = "all",
                         distance = "distance",
                         weights = "weights",
                         subclass = "subclass",
                         data = NULL,
                         include.s.weights = TRUE,
                         drop.unmatched = TRUE)
  print(table(duplicated(matchdata$n_eid)))
  print(table(matchdata$out)) # 匹配后
  print(table(data$out))      # 匹配前
  
  a=table(matchdata$subclass);print(table(a))
  
  print(dim(matchdata))
  outpath = paste0(paste0(paste0('./matched/matched.data/',exp),paste0('.',outcomes[outi])),'.txt')
  write.table(matchdata,outpath,row.names = FALSE,sep=' ',quote = FALSE)
}

## 开始估计smoothing term P-value

var_cal_binomial = function(object1,object2,y,i) # A用glm算，Wei用o算
  # object1: stage 1 model
  # object2: gam model
  # i: No. of the interested smoothing term
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
  
  #Wei = diag(object2$weights)
  S = (object2$smooth[[1]])$S[[1]]
  S = rbind(matrix(rep(0,ncol(Xp)*(ncol(Xp)-nrow(S))),nrow = (ncol(Xp)-nrow(S))),
            cbind(matrix(rep(0,(ncol(Xp)-nrow(S))*nrow(S)),ncol = (ncol(Xp)-nrow(S))),S))
  lambda = object2$sp[[1]]
  W = Xp
  #V1 = solve(t(W)%*%Wei%*%W+lambda*S) %*% crossprod(A) %*% solve(t(W)%*%Wei%*%W)
  V = vcov(object2) %*% crossprod(A) %*% solve(solve(vcov(object2))-lambda*S) # 这样算快一些
  
  return(V)
}

var_cal_binomial = function(object1,object2,y,i) 
  # object1: stage 1 model
  # object2: gam model
  # i: No. of the interested smoothing term
{
  Xp = model.matrix(object2)
  #c2 = glm(y ~ Xp-1, family = "binomial",maxit=50)
  B = matrix(coef(object2),ncol=1)
  etahat = object2$linear.predictors
  miuhat = object2$fitted.values
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
  
  #Wei = diag(object2$weights)
  S = (object2$smooth[[1]])$S[[1]]
  S = rbind(matrix(rep(0,ncol(Xp)*(ncol(Xp)-nrow(S))),nrow = (ncol(Xp)-nrow(S))),
            cbind(matrix(rep(0,(ncol(Xp)-nrow(S))*nrow(S)),ncol = (ncol(Xp)-nrow(S))),S))
  lambda = object2$sp[[1]]
  W = Xp
  V = vcov(object2) %*% crossprod(A) %*% solve(solve(vcov(object2))-lambda*S) 
  
  return(V)
}


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



explist = c('BMI','WC','HIP','WHR','BF','BMR','DBP','SBP','sleep_duration','alcohol','coffee')
outcomes = c('I11','I20','I21','I70','I251','I48','I50','I61','I63')
p = data.frame()

for (expi in 1:length(explist)) {
  exp = explist[expi]
  print(exp)
  for (outi in 1:length(outcomes)) {
    print(outi)
    input = paste0(paste0(paste0('./matched/matched.data/',exp),paste0('.',outcomes[outi])),'.txt')
    data = as.data.frame(fread(input))
    data$exp = data[,which(colnames(data)==exp)]
    c1 = lm(exp ~ PRS
            +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data)
    r = c1$residuals
    data$r = r
    
    data$out = data[,which(colnames(data)==outcomes[outi])]
    o = gam(out ~ s(exp,bs="cr")+r+
              pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data,family = binomial())
    
    V = var_cal_binomial(c1,o,data$out,1)
    p[expi,outi] = p_cal(o,V,1)
  }
}
write.csv(p,'./matched/P_corrected.csv',row.names = FALSE,quote = FALSE)



## 画图
input = './matched/matched.data/coffee.'

data = fread(paste0(input,'I11.txt')) 
c1 = lm(coffee ~ PRS+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data)
r = c1$residuals; data$r = r
c_I11 = gam(I11 ~ s(coffee,bs="cr")+r+
              pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data,family = binomial())


data = fread(paste0(input,'I20.txt'))
c1 = lm(coffee ~ PRS
        +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data)
r = c1$residuals; data$r = r
c_I20 = gam(I20 ~ s(coffee,bs="cr")+r+
              pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data,family = binomial())


data = fread(paste0(input,'I21.txt'))
c1 = lm(coffee ~ PRS
        +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data)
r = c1$residuals; data$r = r
c_I21 = gam(I21 ~ s(coffee,bs="cr")+r+
              pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data,family = binomial())


data = fread(paste0(input,'I70.txt'))
c1 = lm(coffee ~ PRS
        +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data)
r = c1$residuals; data$r = r
c_I70 = gam(I70 ~ s(coffee,bs="cr")+r+
              pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data,family = binomial())

data = fread(paste0(input,'I251.txt'))
c1 = lm(coffee ~ PRS
        +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data)
r = c1$residuals; data$r = r
c_I251 = gam(I251 ~ s(coffee,bs="cr")+r+
               pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data,family = binomial())

data = fread(paste0(input,'I48.txt'))
c1 = lm(coffee ~ PRS
        +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data)
r = c1$residuals; data$r = r
c_I48 = gam(I48 ~ s(coffee,bs="cr")+r+
              pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data,family = binomial())

data = fread(paste0(input,'I50.txt'))
c1 = lm(coffee ~ PRS
        +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data)
r = c1$residuals; data$r = r
c_I50 = gam(I50 ~ s(coffee,bs="cr")+r+
              pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data,family = binomial())

data = fread(paste0(input,'I61.txt'))
c1 = lm(coffee ~ PRS
        +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data)
r = c1$residuals; data$r = r
c_I61 = gam(I61 ~ s(coffee,bs="cr")+r+
              pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data,family = binomial())

data = fread(paste0(input,'I63.txt'))
c1 = lm(coffee ~ PRS
        +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data)
r = c1$residuals; data$r = r
c_I63 = gam(I63 ~ s(coffee,bs="cr")+r+
              pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data,family = binomial())






par(mfrow=c(3,3))
plot.gam(c_I11,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         family = "serif",cex.lab = 1.3,cex.axis=1.3)
title(main = expression('Hypertensive heart disease'), family = "serif", cex.main = 1.3, font.main= 5)

plot.gam(c_I20,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         family = "serif",cex.lab = 1.3,cex.axis=1.3)
title(main = expression('Angina pectoris'), family = "serif", cex.main = 1.3, font.main= 5)

plot.gam(c_I21,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         family = "serif",cex.lab = 1.3,cex.axis=1.3)
title(main = expression('Acute myocardial infarction'), family = "serif", cex.main = 1.3, font.main= 5)

plot.gam(c_I70,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         family = "serif",cex.lab = 1.3,cex.axis=1.3)
title(main = expression('Atherosclerosis'), family = "serif", cex.main = 1.3, font.main= 5)

plot.gam(c_I251,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         family = "serif",cex.lab = 1.3,cex.axis=1.3)
title(main = expression('Atherosclerotic heart disease'), family = "serif", cex.main = 1.3, font.main= 5)

plot.gam(c_I48,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         family = "serif",cex.lab = 1.3,cex.axis=1.3)
title(main = expression('Atrial fibrillation and flutter'), family = "serif", cex.main = 1.3, font.main= 5)

plot.gam(c_I50,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         family = "serif",cex.lab = 1.3,cex.axis=1.3)
title(main = expression('Heart failure'), family = "serif", cex.main = 1.3, font.main= 5)

plot.gam(c_I61,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         family = "serif",cex.lab = 1.3,cex.axis=1.3)
title(main = expression('Intracerebral haemorrhage'), family = "serif", cex.main = 1.3, font.main= 5)

plot.gam(c_I63,rug = TRUE,se=TRUE,n=100,jit = TRUE,shade = TRUE,shade.col = "grey90", 
         family = "serif",cex.lab = 1.3,cex.axis=1.3)
title(main = expression('Cerebral infarction'), family = "serif", cex.main = 1.3, font.main= 5)



##  统计一些基线信息
explist = c('BMI','WC','HIP','WHR','BF','BMR','DBP','SBP','sleep_duration','alcohol','coffee')
outcomes = c('I11','I20','I21','I70','I251','I48','I50','I61','I63')
p = data.frame()

for (expi in 1:length(explist)) {
  exp = explist[expi]
  for (outi in 1:length(outcomes)) {
    input = paste0(paste0(paste0('./matched/matched.data/',exp),paste0('.',outcomes[outi])),'.txt')
    data = as.data.frame(fread(input))
    
    p[expi,outi] = table(data$sex)[[1]] / (table(data$sex)[[1]]+table(data$sex)[[2]])
    
  }
}
write.csv(p,'./matched/samplesize.csv',row.names = FALSE,quote = FALSE)


## 工具变量所解释的暴露方差比例
explist = c('BMI','WC','HIP','WHR','BF','BMR','DBP','SBP','sleep_duration','alcohol','coffee')
outcomes = c('I11','I20','I21','I70','I251','I48','I50','I61','I63')
r2 = data.frame()
F.value = data.frame()
library(TwoSampleMR)
for (expi in 1:length(explist)) {
  exp = explist[expi]
  for (outi in 1:length(outcomes)) {
    input = paste0(paste0(paste0('./matched/matched.data/',exp),paste0('.',outcomes[outi])),'.txt')
    data = as.data.frame(fread(input))
    
    data$exp = data[,which(colnames(data)==exp)]
    c1 = lm(exp ~ PRS
            +pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10,data = data)
    
    r2[expi,outi] = summary(c1)$"r.squared"
    F.value[expi,outi] = summary(c1)$fstatistic[[1]]
  }
}
write.csv(r2,'./matched/r2.csv',row.names = FALSE,quote = FALSE)
write.csv(F.value,'./matched/F.value.csv',row.names = FALSE,quote = FALSE)


## 线性MR
explist = c('BMI','WC','HIP','WHR','BF','BMR','DBP','SBP','sleep_duration','alcohol','coffee')
outcomes = c('I11','I20','I21','I70','I251','I48','I50','I61','I63')
p_sp = data.frame()
p_sr = data.frame()
library(OneSampleMR)
for (expi in 1:length(explist)) {
  exp = explist[expi]
  for (outi in 1:length(outcomes)) {
    input = paste0(paste0(paste0('./matched/matched.data/',exp),paste0('.',outcomes[outi])),'.txt')
    data = as.data.frame(fread(input))
    
    data$exp = data[,which(colnames(data)==exp)]
    data$out = data[,which(colnames(data)==outcomes[outi])]
    
    sp = tsps(out ~ exp | PRS,data = data,link = "logit")
    l = sp$estci[4,2];u = sp$estci[4,3];b = sp$estci[4,1];se = (u-b)/1.96
    p_sp[expi,outi] = pnorm(-abs(b/se))*2
    
    sr = tsri(out ~ exp | PRS,data = data,link = "logit")
    l = sr$estci[4,2];u = sr$estci[4,3];b = sr$estci[4,1];se = (u-b)/1.96
    p_sr[expi,outi] = pnorm(-abs(b/se))*2
  }
}
write.csv(p_sp,'./matched/p_sp.csv',row.names = FALSE,quote = FALSE)







### Over ###

