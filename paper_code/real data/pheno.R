

library(data.table)
data = fread('data1.csv')
names(data)

icd10 = data[,1:244]
icd9 = data[,c(1,245:291)]


BMI = data[,c('n_eid','n_21001_0_0')]
WC = data[,c('n_eid','n_48_0_0')]
HIP = data[,c('n_eid','n_49_0_0')]
BF = data[,c('n_eid','n_23099_0_0')]
BMR = data[,c('n_eid','n_23105_0_0')]
alcohol = data[,c("n_eid","n_1558_0_0","n_20117_0_0","n_4407_0_0","n_4418_0_0","n_4429_0_0","n_4440_0_0","n_4451_0_0","n_4462_0_0",
                  "n_1568_0_0","n_1578_0_0","n_1588_0_0","n_1598_0_0","n_1608_0_0","n_5364_0_0")]
coffee = data[,c("n_eid","n_1498_0_0","n_1508_0_0")]

names(BMI)[2] = 'BMI'
names(WC)[2] = 'WC'
names(HIP)[2] = 'HIP'
names(BF)[2] = 'BF'
names(BMR)[2] = 'BMR'
names(coffee)[2:3] = c('intake','type')
names(alcohol)[2:15] = c('freq','status',
                         'mon_red','mon_cham_white','mon_beer_cider','mon_spirits','mon_fortified','mon_other',
                         'week_red','week_cham_white','week_beer_cider','week_spirits','week_fortified','week_other')

WHR = data[,c('n_eid','n_48_0_0','n_49_0_0')]
WHR$WHR = WHR$n_48_0_0/WHR$n_49_0_0
WHR = WHR[,c('n_eid','WHR')]

table(alcohol$freq)
alcohol = alcohol[-which(alcohol$freq==-3 | alcohol$freq==6),] # 461181 

week_alcohol = alcohol[which(alcohol$freq<4),c('n_eid','freq','status',
                                               'week_red','week_cham_white','week_beer_cider','week_spirits','week_fortified')]
table(week_alcohol$freq) # 346446
table(week_alcohol$status)
table(week_alcohol[,3])
week_alcohol = week_alcohol[-which(week_alcohol$week_red<0 | week_alcohol$week_cham_white<0 | week_alcohol$week_beer_cider<0 | 
                                     week_alcohol$week_spirits<0 | week_alcohol$week_fortified<0),] # 340230

week_alcohol$week_amount = week_alcohol$week_red * 2 + 
                      week_alcohol$week_cham_white * 2 + 
                      week_alcohol$week_beer_cider * 2 +
                      week_alcohol$week_spirits * 1 +
                      week_alcohol$week_fortified * 2
summary(week_alcohol$week_amount)

mon_alcohol = alcohol[which(alcohol$freq==4 | alcohol$freq==5),c('n_eid','freq','status',
                                                                 'mon_red','mon_cham_white','mon_beer_cider','mon_spirits','mon_fortified','mon_other')]
# 113836
table(mon_alcohol$freq)
table(mon_alcohol$status)
table(mon_alcohol[,5])

mon_alcohol = mon_alcohol[-which(mon_alcohol$mon_red<0 | mon_alcohol$mon_cham_white<0 | mon_alcohol$mon_beer_cider<0 | 
                                   mon_alcohol$mon_spirits<0 | mon_alcohol$mon_fortified<0),] # 112613

mon_alcohol = mon_alcohol[-which(is.na(mon_alcohol$mon_red) | is.na(mon_alcohol$mon_cham_white) | 
                                   is.na(mon_alcohol$mon_beer_cider) | 
                                   is.na(mon_alcohol$mon_spirits) | is.na(mon_alcohol$mon_fortified)),] # 39595

a = mon_alcohol[which(is.na(mon_alcohol$mon_red) | is.na(mon_alcohol$mon_cham_white) | 
                                   is.na(mon_alcohol$mon_beer_cider) | 
                                   is.na(mon_alcohol$mon_spirits) | is.na(mon_alcohol$mon_fortified)),]  # 去掉的这一部分这几个变量全是NA



mon_alcohol$mon_amount = mon_alcohol$mon_red * 2 + 
  mon_alcohol$mon_cham_white * 2 + 
  mon_alcohol$mon_beer_cider * 2 +
  mon_alcohol$mon_spirits * 1 +
  mon_alcohol$mon_fortified * 2
mon_alcohol$week_amount = mon_alcohol$mon_amount/4.3
head(week_alcohol)
head(mon_alcohol)
a = week_alcohol[,c('n_eid','week_amount')]
b = mon_alcohol[,c('n_eid','week_amount')]
alcohol_new = rbind(a,b) #379825

write.table(BMI,'BMI.txt',sep=' ',row.names = FALSE,quote = FALSE)
write.table(WC,'WC.txt',sep=' ',row.names = FALSE,quote = FALSE)
write.table(HIP,'HIP.txt',sep=' ',row.names = FALSE,quote = FALSE)
write.table(BF,'BF.txt',sep=' ',row.names = FALSE,quote = FALSE)
write.table(BMR,'BMR.txt',sep=' ',row.names = FALSE,quote = FALSE)
write.table(coffee,'coffee.txt',sep=' ',row.names = FALSE,quote = FALSE)
write.table(alcohol_new,'alcohol.txt',sep=' ',row.names = FALSE,quote = FALSE)
write.table(WHR,'WHR.txt',sep=' ',row.names = FALSE,quote = FALSE)




self = fread('self-reported.csv')
icd10 = data[,1:244]
icd10 = data.frame(icd10)
for (i in 2:ncol(icd10)) {
  icd10[which(icd10[,i]=='I11'),'I11'] = 1
  icd10[which(icd10[,i]=='I110'),'I11'] = 1 # 130
  icd10[which(icd10[,i]=='I119'),'I11'] = 1 # 279
  
  
  icd10[which(icd10[,i]=='I20'),'I20'] = 1
  icd10[which(icd10[,i]=='I200'),'I20'] = 1 # 8195
  icd10[which(icd10[,i]=='I201'),'I20'] = 1 # 352
  icd10[which(icd10[,i]=='I208'),'I20'] = 1 # 3066
  icd10[which(icd10[,i]=='I209'),'I20'] = 1 # 28112
  
  icd10[which(icd10[,i]=='I21'),'I21'] = 1
  icd10[which(icd10[,i]=='I210'),'I21'] = 1 # 2953
  icd10[which(icd10[,i]=='I211'),'I21'] = 1 # 3588
  icd10[which(icd10[,i]=='I212'),'I21'] = 1 # 479
  icd10[which(icd10[,i]=='I213'),'I21'] = 1  # 343
  icd10[which(icd10[,i]=='I214'),'I21'] = 1 # 5261
  icd10[which(icd10[,i]=='I219'),'I21'] = 1 # 5461
  icd10[which(icd10[,i]=='I21X'),'I21'] = 1 # 0
  
  icd10[which(icd10[,i]=='I70'),'I70'] = 1  # 0
  icd10[which(icd10[,i]=='I700'),'I70'] = 1 # 299
  icd10[which(icd10[,i]=='I7000'),'I70'] = 1  # 959
  icd10[which(icd10[,i]=='I7001'),'I70'] = 1 # 13
  icd10[which(icd10[,i]=='I701'),'I70'] = 1 # 120
  icd10[which(icd10[,i]=='I7010'),'I70'] = 1  # 157
  icd10[which(icd10[,i]=='I7011'),'I70'] = 1 # 4
  icd10[which(icd10[,i]=='I702'),'I70'] = 1 # 675
  icd10[which(icd10[,i]=='I7020'),'I70'] = 1  # 1267
  icd10[which(icd10[,i]=='I7021'),'I70'] = 1 # 238
  icd10[which(icd10[,i]=='I708'),'I70'] = 1  # 91
  icd10[which(icd10[,i]=='I7080'),'I70'] = 1  # 265
  icd10[which(icd10[,i]=='I7081'),'I70'] = 1  # 14
  icd10[which(icd10[,i]=='I709'),'I70'] = 1  # 106
  icd10[which(icd10[,i]=='I7090'),'I70'] = 1  # 122
  icd10[which(icd10[,i]=='I7091'),'I70'] = 1  # 2
  
  icd10[which(icd10[,i]=='I251'),'I251'] = 1 # 35532
  
  icd10[which(icd10[,i]=='I48'),'I48'] = 1  # 21778
  icd10[which(icd10[,i]=='I480'),'I48'] = 1 # 6773
  icd10[which(icd10[,i]=='I481'),'I48'] = 1  # 1315
  icd10[which(icd10[,i]=='I482'),'I48'] = 1  # 595
  icd10[which(icd10[,i]=='I483'),'I48'] = 1  # 220
  icd10[which(icd10[,i]=='I484'),'I48'] = 1  # 92
  icd10[which(icd10[,i]=='I489'),'I48'] = 1  # 22529
  
  icd10[which(icd10[,i]=='I50'),'I50'] = 1
  icd10[which(icd10[,i]=='I500'),'I50'] = 1 # 7407
  icd10[which(icd10[,i]=='I501'),'I50'] = 1  # 9144
  icd10[which(icd10[,i]=='I509'),'I50'] = 1  # 6582
  
  icd10[which(icd10[,i]=='I61'),'I61'] = 1 
  icd10[which(icd10[,i]=='I610'),'I61'] = 1 # 208
  icd10[which(icd10[,i]=='I611'),'I61'] = 1 # 333
  icd10[which(icd10[,i]=='I612'),'I61'] = 1  # 109
  icd10[which(icd10[,i]=='I613'),'I61'] = 1  # 60
  icd10[which(icd10[,i]=='I614'),'I61'] = 1  # 114
  icd10[which(icd10[,i]=='I615'),'I61'] = 1 # 210
  icd10[which(icd10[,i]=='I616'),'I61'] = 1  # 61
  icd10[which(icd10[,i]=='I618'),'I61'] = 1  # 205
  icd10[which(icd10[,i]=='I619'),'I61'] = 1  # 989
  
  icd10[which(icd10[,i]=='I63'),'I63'] = 1
  icd10[which(icd10[,i]=='I630'),'I63'] = 1 # 85
  icd10[which(icd10[,i]=='I631'),'I63'] = 1 # 22
  icd10[which(icd10[,i]=='I632'),'I63'] = 1 # 305
  icd10[which(icd10[,i]=='I633'),'I63'] = 1  # 710
  icd10[which(icd10[,i]=='I634'),'I63'] = 1  # 444
  icd10[which(icd10[,i]=='I635'),'I63'] = 1 # 1937
  icd10[which(icd10[,i]=='I636'),'I63'] = 1 # 13
  icd10[which(icd10[,i]=='I638'),'I63'] = 1  # 603
  icd10[which(icd10[,i]=='I639'),'I63'] = 1  # 5794

}
table(icd10$I11) # 397
table(icd10$I20) # 32256
table(icd10$I21) # 16283
table(icd10$I70) # 3567
table(icd10$I251)  # 35532
table(icd10$I48)  # 36321
table(icd10$I50) # 16986
table(icd10$I61) # 1956
table(icd10$I63) # 8936