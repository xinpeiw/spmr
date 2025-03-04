
qc = fread('./pheo_data/sample_QC.txt')
names(qc)[2:5] = c('ethnic','sex_g','kinship','outlier_hete_missing')
cov = fread('./pheo_data/covariates_adj.txt')
names(cov)[2:15] = c('sex','age','centre','batch','pc1','pc2','pc3','pc4','pc5','pc6','pc7','pc8','pc9','pc10')

## QC ##
sex_s = cov[,c('n_eid','sex')]
qc = merge(qc,sex_s,by='n_eid')

wB = qc[which(qc$ethnic=='1001'),'n_eid'] 
kinship = qc[-which(qc$kinship=='-1' | qc$kinship=='10'),'n_eid']
outlier = qc[-which(qc$outlier_hete_missing=='1'),'n_eid']
sex_d = qc[-which(qc$sex_g != qc$sex | is.na(qc$sex_g)),'n_eid']

qc_pass =  intersect(intersect(intersect(wB$n_eid,kinship$n_eid),outlier$n_eid),sex_d$n_eid)  # 429655


cov_age = cov[-which(is.na(cov$age)),] # 502408
cov_batch = cov_age[-which(is.na(cov_age$batch)),] # 488170


explist = c('BMI','WC','HIP','WHR','BF','BMR','DBP','SBP','sleep_duration','alcohol','coffee')
for (i in 1:length(explist)) {
  exp = explist[i]
  
  prspath = paste0(paste0("prs_",exp),'.txt')
  exppath = paste0(paste0("./pheo_data/exposure/",exp),'.txt')
  outcomepath = "./pheo_data/outcome.txt"
  
  prs = fread(prspath)
  exposure = fread(exppath)
  outcome = fread(outcomepath)
  
  exposure = exposure[-which(is.na(exposure[,2])),]
  d = merge(merge(prs,exposure,by='n_eid'),outcome,by='n_eid')
  
  d = d[which(d$n_eid %in% qc_pass),] 
  
  d = d[which(d$n_eid %in% cov_batch$n_eid),]
  
  data = merge(d,cov_batch,by='n_eid')   
  outpath = paste0(paste0("data_",exp),'.txt')
  write.table(data,outpath,sep=' ',row.names = FALSE,quote = FALSE)
  
  print(exp)
  print(dim(data))
}



exp = 'SBP'

prspath = paste0(paste0("prs_",exp),'.txt')
exppath = paste0(paste0("./pheo_data/exposure/",exp),'.txt')
outcomepath = "./pheo_data/outcome.txt"

prs = fread(prspath)
exposure = fread(exppath)
outcome = fread(outcomepath)

exposure = exposure[-which(is.na(exposure$n_4080_0_0) & is.na(exposure$n_4080_0_1)),] 

e1 = exposure[which(is.na(exposure$n_4080_0_0)),]
e1$SBP = e1$n_4080_0_1
e2 = exposure[which(is.na(exposure$n_4080_0_1)),] 
e2$SBP = e2$n_4080_0_0
e3 = exposure[-which(is.na(exposure$n_4080_0_0) | is.na(exposure$n_4080_0_1)),]
e3$SBP = (e3$n_4080_0_0+e3$n_4080_0_1)/2
exposure = rbind(e1,e2,e3)

d = merge(merge(prs,exposure,by='n_eid'),outcome,by='n_eid') 

d = d[which(d$n_eid %in% qc_pass),]  

d = d[which(d$n_eid %in% cov_batch$n_eid),]

data = merge(d,cov_batch,by='n_eid') 

outpath = paste0(paste0("data_",exp),'.txt')
write.table(data,outpath,sep=' ',row.names = FALSE,quote = FALSE)







