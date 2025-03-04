
library(data.table)
library(TwoSampleMR)


iv_select = function(data) {
  b = data[which(data$pval<5e-8),]
  b_format = format_data(
    b,
    type = "exposure",
    snps = NULL,
    header = TRUE,
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    samplesize_col = "samplesize",
    min_pval = 1e-200
  )
  
  b_clump = clump_data(
    b_format,
    clump_kb = 10000,
    clump_r2 = 0.001,
    clump_p1 = 1,
    clump_p2 = 1,
    pop = "EUR"
  )
  
  b_iv = b_clump[,c('SNP','effect_allele.exposure','other_allele.exposure','beta.exposure')]
  
  return(b_iv)
}

## BMI ##
a = fread("BMI_data.txt")
head(a)
names(a)[3:10] = c('SNP','effect_allele','other_allele','eaf','beta','se','pval','samplesize')
a_iv = iv_select(a)
dim(a_iv)
write.table(a_iv,'IV_BMI.txt',row.names = FALSE,sep=' ',quote = FALSE) # 521个

## WC ##
a = fread("GIANT_2015_WC_COMBINED_EUR.txt")
head(a)
names(a) = c('SNP','effect_allele','other_allele','eaf','beta','se','pval','samplesize')
a_iv = iv_select(a)
dim(a_iv)
write.table(a_iv,'IV_WC.txt',row.names = FALSE,sep=' ',quote = FALSE)

## HIP ##
a = fread("GIANT_2015_HIP_COMBINED_EUR.txt")
head(a)
names(a) = c('SNP','effect_allele','other_allele','eaf','beta','se','pval','samplesize')
a_iv = iv_select(a)
dim(a_iv)
write.table(a_iv,'IV_HIP.txt',row.names = FALSE,sep=' ',quote = FALSE)

## WHR ##
a = fread("GIANT_2015_WHR_COMBINED_EUR.txt")
head(a)
a = a[,-c(2,3)]
names(a) = c('SNP','effect_allele','other_allele','eaf','beta','se','pval','samplesize')
a_iv = iv_select(a)
dim(a_iv)
write.table(a_iv,'IV_WHR.txt',row.names = FALSE,sep=' ',quote = FALSE)

## BF ##
a = fread("26833246-GCST003435-EFO_0007800.h.tsv")
head(a)
a = a[,c('variant_id','effect_allele','other_allele','effect_allele_frequency','beta','standard_error','p_value','n')]
names(a) = c('SNP','effect_allele','other_allele','eaf','beta','se','pval','samplesize')
a_iv = iv_select(a)
dim(a_iv)
write.table(a_iv,'IV_BF.txt',row.names = FALSE,sep=' ',quote = FALSE)

## DBP ##
a = fread("Evangelou_30224653_DBP.txt")  
library("tidyr")
b = separate(a,'MarkerName',c('chr','pos','type'),sep = ":")
b = b[which(b$P<5e-8),]
all = fread('variants.tsv')
head(all)
all = all[,c('chr','pos','rsid')]
all$chr = as.numeric(all$chr); all$pos = as.numeric(all$pos)
b$chr = as.numeric(b$chr); b$pos = as.numeric(b$pos)
c = merge(b,all,by=c('chr','pos'))
head(c)
c = c[,c('rsid','Allele1','Allele2','Freq1','Effect','StdErr','P','TotalSampleSize')]
names(c) = c('SNP','effect_allele','other_allele','eaf','beta','se','pval','samplesize')
a_iv = iv_select(c)
dim(a_iv)
write.table(a_iv,'IV_DBP.txt',row.names = FALSE,sep=' ',quote = FALSE)


## SBP ##
a = fread("Evangelou_30224653_SBP.txt")  
library("tidyr")
a = a[which(a$P<5e-8),]
b = separate(a,'MarkerName',c('chr','pos','type'),sep = ":")
all = fread('variants.tsv')
head(all)
all = all[,c('chr','pos','rsid')]
all$chr = as.numeric(all$chr); all$pos = as.numeric(all$pos)
b$chr = as.numeric(b$chr); b$pos = as.numeric(b$pos)
c = merge(b,all,by=c('chr','pos'))
head(c)
c = c[,c('rsid','Allele1','Allele2','Freq1','Effect','StdErr','P','TotalSampleSize')]
names(c) = c('SNP','effect_allele','other_allele','eaf','beta','se','pval','samplesize')
a_iv = iv_select(c)
dim(a_iv)
write.table(a_iv,'IV_SBP.txt',row.names = FALSE,sep=' ',quote = FALSE)


## Sleep duration ##
a = fread("sleepdurationsumstats.txt")
head(a)
# a = a[,c('SNP','ALLELE1','ALLELE0','A1FREQ','BETA_SLEEPDURATION','SE_SLEEPDURATION','P_SLEEPDURATION')]
names(a) = c('SNP','effect_allele','other_allele','eaf','beta','se','pval')
a_iv = iv_select(a)
dim(a_iv)
write.table(a_iv,'IV_sleep_duration.txt',row.names = FALSE,sep=' ',quote = FALSE)


## alcohol ##
a = fread("DrinksPerWeek.txt")
head(a)
a = a[,c('RSID','ALT','REF','AF','BETA','SE','PVALUE','N')]
names(a) = c('SNP','effect_allele','other_allele','eaf','beta','se','pval','samplesize')
a_iv = iv_select(a)
dim(a_iv)
write.table(a_iv,'IV_alcohol.txt',row.names = FALSE,sep=' ',quote = FALSE) 


## coffee ##
a = fread("coffee.assess.nobmi")
head(a)
a = a[,c('ID','ALT','REF','ALT_FREQ','BETA','SE','P')]
names(a) = c('SNP','effect_allele','other_allele','eaf','beta','se','pval')
a_iv = iv_select(a)
dim(a_iv)
write.table(a_iv,'IV_coffee.txt',row.names = FALSE,sep=' ',quote = FALSE) 


## BMR ##
a = fread("29892013-GCST90029025-EFO_0007777-Build37.tsv")
head(a)
a = a[,c('variant_id','effect_allele','other_allele','effect_allele_frequency','beta','standard_error','p_value')]
names(a) = c('SNP','effect_allele','other_allele','eaf','beta','se','pval')
a_iv = iv_select(a)
dim(a_iv)
write.table(a_iv,'./IV_BMR.txt',row.names = FALSE,sep=' ',quote = FALSE) 



folder_path <- getwd()
file_names <- list.files(folder_path)
iv_list = c()
for (i in 1:length(file_names)) {
  dat = fread(file_names[i])
  if (i==1) {
    iv_list = dat$SNP
  } else {
    iv_list = union(iv_list,dat$SNP)
  }
}

length(iv_list)
length(unique(iv_list))
write.table(iv_list,'iv_list.txt',row.names = FALSE,sep = ' ',quote = FALSE)


