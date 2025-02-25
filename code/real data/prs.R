setwd('C:/ShareCache/王鑫培_2111110203/AAA小论文/REAL DATA')

myfile = list.files('./snp')   
rawfile = myfile[grep(myfile, pattern =".raw$")]   
rawfile

# 先形成UKB中snp的rsid,count_allele,noncount_allele
snplist = c()
for (i in 1:length(rawfile)) {
  data = read.table(paste0('./snp/',rawfile[i]),header = TRUE)
  s = colnames(data)[7:ncol(data)]
  snplist = c(snplist,s)
}
length(snplist)
df = data.frame('combine'=NA,'rsID'=NA,'count'=NA,'noncount'=NA)
for (i in 1:length(snplist)) {
  s = strsplit(snplist[i],split = '_')[[1]]
  s2 = strsplit(s[2],split = '\\..')[[1]]
  s3 = strsplit(s2[2],split = '\\.')
  df[i,'combine'] = snplist[i]
  df[i,'rsID'] = s[1]
  df[i,'count'] = s2[1]
  df[i,'noncount'] = s3
}


## 把每个暴露的IV基因型信息提取出来，生成raw文件 ##
## BMI ##
IV_BMI = fread('./GWAS summary data/IV/IV_BMI.txt')
df_bmi = df[which(df$rsID %in% IV_BMI$SNP),]
dim(df_bmi); dim(IV_BMI)

for (i in 1:length(rawfile)) {
  data = read.table(paste0('./snp/',rawfile[i]),header = TRUE)
  a1 = data[,1:6]
  a2 = data[,colnames(data) %in% df_bmi$combine]
  a = cbind(a1,a2)
  if (i==1) {
    raw_bmi = a
  } else {
    raw_bmi = merge(raw_bmi,a,by=c('FID','IID','PAT','MAT','SEX','PHENOTYPE'))
  }
}

dim(raw_bmi); dim(df_bmi)
write.table(raw_bmi,'raw_bmi.txt',row.names = FALSE,sep=' ',quote = FALSE)


## WC ##
IV_WC = fread('./GWAS summary data/IV/IV_WC.txt')
df_WC = df[which(df$rsID %in% IV_WC$SNP),]
dim(df_WC); dim(IV_WC)

for (i in 1:length(rawfile)) {
  data = read.table(paste0('./snp/',rawfile[i]),header = TRUE)
  a1 = data[,1:6]
  in_col = colnames(data)[colnames(data) %in% df_WC$combine]
  a2 = subset(data,select = in_col)
  colnames(a2)
  a = cbind(a1,a2)
  if (i==1) {
    raw_WC = a
  } else {
    raw_WC = merge(raw_WC,a,by=c('FID','IID','PAT','MAT','SEX','PHENOTYPE'))
  }
  colnames(raw_WC)
}

dim(raw_WC); dim(df_WC)
write.table(raw_WC,'raw_WC.txt',row.names = FALSE,sep=' ',quote = FALSE)


## 写个函数干这件事吧
library(data.table)
inte = function(path) {
  IV_data = fread(path)
  df_data = df[which(df$rsID %in% IV_data$SNP),]
  print(dim(df_data)); print(dim(IV_data))
  
  for (i in 1:length(rawfile)) {
    data = read.table(paste0('./snp/',rawfile[i]),header = TRUE)
    a1 = data[,1:6]
    in_col = colnames(data)[colnames(data) %in% df_data$combine]
    a2 = subset(data,select = in_col)
    a = cbind(a1,a2)
    if (i==1) {
      raw_data = a
    } else {
      raw_data = merge(raw_data,a,by=c('FID','IID','PAT','MAT','SEX','PHENOTYPE'))
    }
  }
  print(dim(raw_data))
  return(raw_data)
}

inte('./GWAS summary data/IV/IV_WC.txt')
raw_HIP = inte('./GWAS summary data/IV/IV_HIP.txt')
write.table(raw_HIP,'raw_HIP.txt',row.names = FALSE,sep=' ',quote = FALSE)

raw_WHR = inte('./GWAS summary data/IV/IV_WHR.txt')
write.table(raw_WHR,'raw_WHR.txt',row.names = FALSE,sep=' ',quote = FALSE)

raw_BF = inte('./GWAS summary data/IV/IV_BF.txt')
write.table(raw_BF,'raw_BF.txt',row.names = FALSE,sep=' ',quote = FALSE)

raw_BMR = inte('./GWAS summary data/IV/IV_BMR.txt')
write.table(raw_BMR,'raw_BMR.txt',row.names = FALSE,sep=' ',quote = FALSE)

raw_DBP = inte('./GWAS summary data/IV/IV_DBP.txt')
write.table(raw_DBP,'raw_DBP.txt',row.names = FALSE,sep=' ',quote = FALSE)

raw_SBP = inte('./GWAS summary data/IV/IV_SBP.txt')
write.table(raw_SBP,'raw_SBP.txt',row.names = FALSE,sep=' ',quote = FALSE)

raw_sleep_duration = inte('./GWAS summary data/IV/IV_sleep_duration.txt')
write.table(raw_sleep_duration,'raw_sleep_duration.txt',row.names = FALSE,sep=' ',quote = FALSE)

raw_alcohol = inte('./GWAS summary data/IV/IV_alcohol.txt')
write.table(raw_alcohol,'raw_alcohol.txt',row.names = FALSE,sep=' ',quote = FALSE)

raw_coffee = inte('./GWAS summary data/IV/IV_coffee.txt')
write.table(raw_coffee,'raw_coffee.txt',row.names = FALSE,sep=' ',quote = FALSE)



## 将count allele转化为risk increase allele ##
path = './GWAS summary data/IV/IV_WC.txt'
IV_data = fread(path)
df_data = df[which(df$rsID %in% IV_data$SNP),]

for (i in 1:nrow(IV_data)) {
  b = IV_data$beta.exposure[i]
  a1 = IV_data$effect_allele.exposure[i]; a2 = IV_data$other_allele.exposure[i]
  if (b>0) {
    next
  } else if (b<0) {
    IV_data$effect_allele.exposure[i] = a2
    IV_data$other_allele.exposure[i] = a1
    IV_data$beta.exposure[i] = -b
  } else {
    print(i)
  }
}

names(IV_data)[1] = 'rsID'
d = merge(df_data,IV_data,by = 'rsID')

# 把beta换为对UKB count allele #
inver = c()
for (i in 1:nrow(d)) {
  if ((d[i,'count'] == d[i,'effect_allele.exposure']) & (d[i,'noncount'] == d[i,'other_allele.exposure'])) {
    next
  } else if ((d[i,'count'] == d[i,'other_allele.exposure']) & (d[i,'noncount'] == d[i,'effect_allele.exposure'])){
    inver = union(inver,i)
  } else {
    print(i)
  }
} 


inver_snp = d[inver,]
inver_snp$inte = paste0(paste(paste(inver_snp$rsID,inver_snp$count,sep='_'),inver_snp$noncount,sep='..'),'.')

gene = data.frame(fread('./raw_WC.txt'))
gene_inver = gene[,which(names(gene)%in% inver_snp$combine)]
gene_inver2 = 2-gene_inver  # 把颠倒的转换过来，现在全部都是risk increase allele
gene1 = gene[,-which(names(gene)%in% inver_snp$inte)]
gene_new = cbind(gene1,gene_inver2)

write.table(gene_new,'gene_info_WC.txt',row.names = FALSE,sep = ' ',quote = FALSE)


## 写个函数干这件事 ##
risk_inc = function(path,pathraw) {
  IV_data = fread(path)
  df_data = df[which(df$rsID %in% IV_data$SNP),]
  
  for (i in 1:nrow(IV_data)) {
    b = IV_data$beta.exposure[i]
    a1 = IV_data$effect_allele.exposure[i]; a2 = IV_data$other_allele.exposure[i]
    if (b>0) {
      next
    } else if (b<0) {
      IV_data$effect_allele.exposure[i] = a2
      IV_data$other_allele.exposure[i] = a1
      IV_data$beta.exposure[i] = -b
    } else {
      print(i)
    }
  }
  
  names(IV_data)[1] = 'rsID'
  d = merge(df_data,IV_data,by = 'rsID')
  
  # 把beta换为对UKB count allele #
  inver = c()
  for (i in 1:nrow(d)) {
    if ((d[i,'count'] == d[i,'effect_allele.exposure']) & (d[i,'noncount'] == d[i,'other_allele.exposure'])) {
      next
    } else if ((d[i,'count'] == d[i,'other_allele.exposure']) & (d[i,'noncount'] == d[i,'effect_allele.exposure'])){
      inver = union(inver,i)
    } else {
      print(i)
    }
  } 
  
  inver_snp = d[inver,]

  gene = data.frame(fread(pathraw))
  gene = gene[,-which(names(gene)=='rs62294352_A..C.')]
  gene_inver = gene[,which(names(gene)%in% inver_snp$combine)]
  gene_inver2 = 2-gene_inver  # 把颠倒的转换过来，现在全部都是risk increase allele
  gene1 = gene[,-which(names(gene)%in% inver_snp$combine)]
  gene_new = cbind(gene1,gene_inver2)
  
  return(gene_new)
}

gene_new_WC = risk_inc('./GWAS summary data/IV/IV_WC.txt','./raw_WC.txt')
write.table(gene_new_WC,'gene_info_WC.txt',row.names = FALSE,sep = ' ',quote = FALSE)

gene_new_BMI = risk_inc('./GWAS summary data/IV/IV_BMI.txt','./raw_BMI.txt')
dim(gene_new_BMI)
write.table(gene_new_BMI,'gene_info_BMI.txt',row.names = FALSE,sep = ' ',quote = FALSE)

gene_new_HIP = risk_inc('./GWAS summary data/IV/IV_HIP.txt','./raw_HIP.txt')
dim(gene_new_HIP)
write.table(gene_new_HIP,'gene_info_HIP.txt',row.names = FALSE,sep = ' ',quote = FALSE)

gene_new_WHR = risk_inc('./GWAS summary data/IV/IV_WHR.txt','./raw_WHR.txt')
dim(gene_new_WHR)
write.table(gene_new_WHR,'gene_info_WHR.txt',row.names = FALSE,sep = ' ',quote = FALSE)

gene_new_BF = risk_inc('./GWAS summary data/IV/IV_BF.txt','./raw_BF.txt')
dim(gene_new_BF)
write.table(gene_new_BF,'gene_info_BF.txt',row.names = FALSE,sep = ' ',quote = FALSE)

gene_new_BMR = risk_inc('./GWAS summary data/IV/IV_BMR.txt','./raw_BMR.txt')
gene_new_BMR = gene_new  # 单独分析了，因为有一个allele不匹配的，删去了
dim(gene_new_BMR)
write.table(gene_new_BMR,'gene_info_BMR.txt',row.names = FALSE,sep = ' ',quote = FALSE)

gene_new_DBP = risk_inc('./GWAS summary data/IV/IV_DBP.txt','./raw_DBP.txt')
gene_new_DBP = gene_new
dim(gene_new_DBP) # 单独分析了，因为有一个allele不匹配的，删去了
write.table(gene_new_DBP,'gene_info_DBP.txt',row.names = FALSE,sep = ' ',quote = FALSE)

gene_new_SBP = risk_inc('./GWAS summary data/IV/IV_SBP.txt','./raw_SBP.txt')
dim(gene_new_SBP)
write.table(gene_new_SBP,'gene_info_SBP.txt',row.names = FALSE,sep = ' ',quote = FALSE)

gene_new_sleep_duration = risk_inc('./GWAS summary data/IV/IV_sleep_duration.txt','./raw_sleep_duration.txt')
dim(gene_new_sleep_duration)
write.table(gene_new_sleep_duration,'gene_info_sleep_duration.txt',row.names = FALSE,sep = ' ',quote = FALSE)

gene_new_alcohol = risk_inc('./GWAS summary data/IV/IV_alcohol.txt','./raw_alcohol.txt')
dim(gene_new_alcohol)
write.table(gene_new_alcohol,'gene_info_alcohol.txt',row.names = FALSE,sep = ' ',quote = FALSE)

gene_new_coffee = risk_inc('./GWAS summary data/IV/IV_coffee.txt','./raw_coffee.txt')
dim(gene_new_coffee)
write.table(gene_new_coffee,'gene_info_coffee.txt',row.names = FALSE,sep = ' ',quote = FALSE)


## 对缺失基因型进行众数填补，并计算PRS
explist = c('BMI','WC','HIP','WHR','BF','BMR','DBP','SBP','sleep_duration','alcohol','coffee')

library(DescTools)  
for (j in 1:length(explist)) {
  inpath = paste0(paste0('gene_info_',explist[j]),'.txt')
  outpath = paste0(paste0('prs_',explist[j]),'.txt')
  
  gene = data.frame(fread(inpath))
  for (i in 7:ncol(gene)) {
    gene[,i][is.na(gene[,i])] = Mode(gene[,i],na.rm=TRUE)[1] 
  }
  #计算PRS
  gene$PRS = apply(gene[,7:ncol(gene)],1,sum)
  print(explist[j])
  print(summary(gene$PRS))
  gene = gene[,c('IID','SEX','PRS')]
  names(gene)[1:2] = c('n_eid','sex_gene')
  write.table(gene,outpath,row.names = FALSE,sep = ' ',quote = FALSE)
  rm(gene)
}




















