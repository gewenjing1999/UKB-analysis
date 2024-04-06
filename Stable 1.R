#Stable1
#基线特征case/control

#加载数据
load("imputed_final_20240120_UKB_database_for_analysis_202801_without_prevalent_cancer_male.rdata")
load("imputed_final_20240120_UKB_database_for_analysis_239598_without_prevalent_cancer_female.rdata")

#加载R包
library(data.table)
library(gmodels) #“gmodels”包的CrossTable()函数


####随访时间####
male$diffyear <- as.numeric(as.Date(male$Censoring_date.y)-as.Date(male$Date_of_attending_assessment_centre.y))/365.25  
summary(male$diffyear)

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#9.413  10.398  11.086  11.153  11.775  13.947 

female$diffyear <- as.numeric(as.Date(female$Censoring_date.y)-as.Date(female$Date_of_attending_assessment_centre.y))/365.25  
summary(female$diffyear)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#9.413  10.407  11.088  11.155  11.767  13.966 

data <- rbind(male,female)
summary(data$diffyear)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#9.413  10.401  11.088  11.154  11.770  13.966 
#During a median follow-up period of 11.09 years (IQR: 10.40–11.77)



#家族史数据整理
#对父亲、母亲、兄弟姐妹的家族史情况进行分别统计
##father字段：20107
##mother字段：20110
##sibling字段：20111
##在ukb中，定义13-前列腺癌，5-乳腺癌，4-肠癌，3-肺癌##
bd <- read.table("family_history.tab", header=TRUE)

#################FH_Prostate#################
##father_FH##
bd$FH_Prostate_fa <- ifelse(apply(bd[,2:41],1,function(x,na.rm=TRUE){sum(grepl('^13$',x))>0}),1,0)
##sibling_FH##
bd$FH_Prostate_sib <- ifelse(apply(bd[,86:133],1,function(x,na.rm=TRUE){sum(grepl('^13$',x))>0}),1,0)
bd$FH_Prostate_total=bd$FH_Prostate_fa+bd$FH_Prostate_sib
table(bd$FH_Prostate_total)
#462324  39194    989

###################FH_Breast###################
##mother_FH##
bd$FH_Breast_mo <- ifelse(apply(bd[,42:85],1,function(x,na.rm=TRUE){sum(grepl('^5$',x))>0}),1,0)
##sibling_FH##
bd$FH_Breast_sib <- ifelse(apply(bd[,86:133],1,function(x,na.rm=TRUE){sum(grepl('^5$',x))>0}),1,0)
bd$FH_Breast_total=bd$FH_Breast_mo+bd$FH_Breast_sib
table(bd$FH_Breast_total)
#448693  51381   2433

#####################FH_Colorectal#####################
##father_FH##
bd$FH_Bowel_fa <- ifelse(apply(bd[,2:41],1,function(x,na.rm=TRUE){sum(grepl('^4$',x))>0}),1,0)
##mother_FH##
bd$FH_Bowel_mo <- ifelse(apply(bd[,42:85],1,function(x,na.rm=TRUE){sum(grepl('^4$',x))>0}),1,0)
##sibling_FH##
bd$FH_Bowel_sib <- ifelse(apply(bd[,86:133],1,function(x,na.rm=TRUE){sum(grepl('^4$',x))>0}),1,0)
bd$FH_Bowel_total=bd$FH_Bowel_fa+bd$FH_Bowel_mo+bd$FH_Bowel_sib
table(bd$FH_Bowel_total)
#446566  52648   3193    100

#####################FH_Lung#####################
##father_FH##
bd$FH_Lung_fa <- ifelse(apply(bd[,2:41],1,function(x,na.rm=TRUE){sum(grepl('^3$',x))>0}),1,0)
##mother_FH##
bd$FH_Lung_mo <- ifelse(apply(bd[,42:85],1,function(x,na.rm=TRUE){sum(grepl('^3$',x))>0}),1,0)
##sibling_FH##
bd$FH_Lung_sib <- ifelse(apply(bd[,86:133],1,function(x,na.rm=TRUE){sum(grepl('^3$',x))>0}),1,0)
bd$FH_Lung_total=bd$FH_Lung_fa+bd$FH_Lung_mo+bd$FH_Lung_sib
table(bd$FH_Lung_total)
#439100  59143   4087    177


re<-bd[,c(1,136,139,143,147)] ####取出刚刚得到的特定肿瘤家族史变量
names(re)
#[1] "f.eid"             "FH_Prostate_total" "FH_Breast_total"   "FH_Bowel_total"    "FH_Lung_total" 

male <- merge(male,re,by="f.eid")
female <- merge(female,re,by="f.eid")

####男性总体肿瘤家族史包含肠癌，肺癌，前列腺癌   FH_total_2男性总体肿瘤家族史二分类
male$FH_total <- male$FH_Bowel_total.y+male$FH_Lung_total.y+male$FH_Prostate_total.y
table(male$FH_total)
#0      1      2      3      4      5      6      7 
#146524  47087   8026   1011    129     22      1      1 

male$FH_total_2 <- ifelse(male$FH_total==0,0,1)
table(male$FH_total_2)
#146524  56277


#男性总体肿瘤家族史（肠癌，肺癌，前列腺癌）四分类  ，FH_total_4男性总体肿瘤家族史四分类 202801
male$FH_total_4=0
male[male$FH_total==1,]$FH_total_4=1
male[male$FH_total==2,]$FH_total_4=2
male[male$FH_total>2,]$FH_total_4=3
table(male$FH_total_4)
#0      1      2      3 
#146524  47087   8026   1164 
#72.25   23.22   3.96   0.57
#有家族史人数：56277
#47087   8026   1164
#83.67   14.26  2.07

#男性肠癌家族史三分类0，1，≥2   FH_Bowel_3
male$FH_Bowel_3=0
male[male$FH_Bowel_total.y==1,]$FH_Bowel_3=1
male[male$FH_Bowel_total.y>1,]$FH_Bowel_3=2
table(male$FH_Bowel_3)
#180048  21386   1367

#男性肺癌家族史三分类0，1，≥2   FH_Lung_3
male$FH_Lung_3=0
male[male$FH_Lung_total.y==1,]$FH_Lung_3=1
male[male$FH_Lung_total.y>1,]$FH_Lung_3=2
table(male$FH_Lung_3)
#178146  22912   1743

#男性前列腺癌家族史三分类0，1，≥2
male$FH_Prostate_3=0
male[male$FH_Prostate_total.y==1,]$FH_Prostate_3=1
male[male$FH_Prostate_total.y>1,]$FH_Prostate_3=2
table(male$FH_Prostate_3)
#187028  15378    395




####女性肿瘤家族史包含肠癌，肺癌，乳腺癌          FH_total_2女性总体肿瘤家族史二分类
female$FH_total <- female$FH_Bowel_total.y+female$FH_Lung_total.y+female$FH_Breast_total.y
table(female$FH_total)
#0      1      2      3      4      5      6      7 
#167534  59088  11281   1498    171     21      4      1 

female$FH_total_2 <- ifelse(female$FH_total==0,0,1)
table(female$FH_total_2)
#167534  72064

#女性总体肿瘤家族史（肠癌，肺癌，乳腺癌）四分类  ，FH_total_4女性总体肿瘤家族史四分类 239598
female$FH_total_4=0
female[female$FH_total==1,]$FH_total_4=1
female[female$FH_total==2,]$FH_total_4=2
female[female$FH_total>2,]$FH_total_4=3
table(female$FH_total_4)
#0      1      2      3 
#167534  59088  11281   1695 
#69.92   24.66  4.71    0.71
#有家族史人数72064
#59088  11281   1695
#81.99  15.65   2.36

#女性肠癌家族史三分类0，1，≥2   FH_Bowel_3
female$FH_Bowel_3=0
female[female$FH_Bowel_total.y==1,]$FH_Bowel_3=1
female[female$FH_Bowel_total.y>1,]$FH_Bowel_3=2
table(female$FH_Bowel_3)
#213800  24426   1372 

#女性肺癌家族史三分类0，1，≥2   FH_Lung_3
female$FH_Lung_3=0
female[female$FH_Lung_total.y==1,]$FH_Lung_3=1
female[female$FH_Lung_total.y>1,]$FH_Lung_3=2
table(female$FH_Lung_3)
#209249  28402   1947


#女性乳腺癌家族史三分类0，1，≥2
female$FH_Breast_3=0
female[female$FH_Breast_total.y==1,]$FH_Breast_3=1
female[female$FH_Breast_total.y>1,]$FH_Breast_3=2
table(female$FH_Breast_3)
#213373  25063   1162


####更新全肿瘤cox生存时间数survival time

#男性 
for (i in 1:nrow(male)) {
  a=min(male[i,"CRC_difftime_new"],male[i,"LC_difftime_new"],male[i,"PRC_difftime_new"],male[i,"BC_difftime_new"],male[i,"GC_difftime_new"],male[i,"PAC_difftime_new"],male[i,"TEC_difftime_new"],
        male[i,"OEC_difftime_new"],male[i,"CEUC_difftime_new"],male[i,"Hodgkin_difftime_new"],male[i,"MM_difftime_new"],male[i,"LOCPC_difftime_new"],male[i,"SKM_difftime_new"],
        male[i,"COUC_difftime_new"],male[i,"BLC_difftime_new"],male[i,"THC_difftime_new"],male[i,"BCNC_difftime_new"],male[i,"KIC_difftime_new"],male[i,"LYLE_difftime_new"],male[i,"OVC_difftime_new"])
  male[i,"survival_time"] <- a
}

#女性 
for (i in 1:nrow(female)) {
  a=min(female[i,"CRC_difftime_new"],female[i,"LC_difftime_new"],female[i,"PRC_difftime_new"],female[i,"BC_difftime_new"],female[i,"GC_difftime_new"],female[i,"PAC_difftime_new"],female[i,"TEC_difftime_new"],
        female[i,"OEC_difftime_new"],female[i,"CEUC_difftime_new"],female[i,"Hodgkin_difftime_new"],female[i,"MM_difftime_new"],female[i,"LOCPC_difftime_new"],female[i,"SKM_difftime_new"],
        female[i,"COUC_difftime_new"],female[i,"BLC_difftime_new"],female[i,"THC_difftime_new"],female[i,"BCNC_difftime_new"],female[i,"KIC_difftime_new"],female[i,"LYLE_difftime_new"],female[i,"OVC_difftime_new"])
  female[i,"survival_time"] <- a
}



data <- rbind(male,female)



#统计经过随访共有多少新发肿瘤病例
table(male$cancer_total_20)
#0      1 
#183133  19668 

table(female$cancer_total_20)
#221575  18023

table(data$cancer_total_20)
#404708  37691


#分别定义case/control
inci=subset(data,data$cancer_total_20==1)
con=subset(data,data$cancer_total_20==0)

male_inci=subset(male,male$cancer_total_20==1)
male_con=subset(male,male$cancer_total_20==0)

female_inci=subset(female,female$cancer_total_20==1)
female_con=subset(female,female$cancer_total_20==0)


###############年龄   # Age was calculated by using the birth date and the date of baseline assessment. #
##male
mean(male_inci$Age)
sd(male_inci$Age)
#[1] 60.87975
#> sd(male_inci$Age)
#[1] 6.309208

mean(male_con$Age)
sd(male_con$Age)
#[1] 55.80693
#> sd(male_con$Age)
#[1] 8.227754

##female
mean(female_inci$Age)
sd(female_inci$Age)
#[1] 58.41092
#> sd(female_inci$Age)
#[1] 7.456502

mean(female_con$Age)
sd(female_con$Age)
#[1] 55.83234
#> sd(female_con$Age)
#[1] 8.030396

##all
mean(inci$Age)
sd(inci$Age)
#[1] 59.69921
#> sd(inci$Age)
#[1] 6.991263

mean(con$Age)
sd(con$Age)
#[1] 55.82084
#> sd(con$Age)
#[1] 8.120296


############height
##male
mean(male_inci$Height_imp)
sd(male_inci$Height_imp)
#> mean(male_inci$Height_imp)
#[1] 175.3077
#> sd(male_inci$Height_imp)
#[1] 6.720223

mean(male_con$Height_imp)
sd(male_con$Height_imp)
#> mean(male_con$Height_imp)
#[1] 175.6694
#> sd(male_con$Height_imp)
#[1] 6.866049


##female
mean(female_inci$Height_imp)
sd(female_inci$Height_imp)
#> mean(female_inci$Height_imp)
#[1] 162.5691
#> sd(female_inci$Height_imp)
#[1] 6.297133

mean(female_con$Height_imp)
sd(female_con$Height_imp)
#> mean(female_con$Height_imp)
#[1] 162.4437
#> sd(female_con$Height_imp)
#[1] 6.314794


##all
mean(inci$Height_imp)
sd(inci$Height_imp)
#> mean(inci$Height_imp)
#[1] 169.2164
#> sd(inci$Height_imp)
#[1] 9.111427

mean(con$Height_imp)
sd(con$Height_imp)
#> mean(con$Height_imp)
#[1] 168.4284
#> sd(con$Height_imp)
#[1] 9.300545


##Normal BMI 
#male
CrossTable(male$BMI_health,male$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  202801 
#| male$cancer_total_20 
#male$BMI_health |         0 |         1 | Row Total | 
#  ----------------|-----------|-----------|-----------|
#  0 |     46903 |      5289 |     52192 | 
# |    1.0965 |   10.2096 |           | 
#  |    0.2561 |    0.2689 |           | 
#  ----------------|-----------|-----------|-----------|
#  1 |    136230 |     14379 |    150609 | 
#  |    0.3800 |    3.5380 |           | 
#  |    0.7439 |    0.7311 |           | 
#  ----------------|-----------|-----------|-----------|
#  Column Total |    183133 |     19668 |    202801 | 
#  |    0.9030 |    0.0970 |           | 
#  ----------------|-----------|-----------|-----------|
# Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  15.22412     d.f. =  1     p =  9.54759e-05 
#Pearson's Chi-squared test with Yates' continuity correction 
#------------------------------------------------------------
#Chi^2 =  15.15723     d.f. =  1     p =  9.891887e-05

#female
CrossTable(female$BMI_health,female$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  239598 
#| female$cancer_total_20 
#female$BMI_health |         0 |         1 | Row Total | 
#  ------------------|-----------|-----------|-----------|
#  0 |     53362 |      4953 |     58315 | 
#  |    5.9496 |   73.1446 |           | 
#  |    0.2408 |    0.2748 |           | 
#  ------------------|-----------|-----------|-----------|
#  1 |    168213 |     13070 |    181283 | 
#  |    1.9139 |   23.5291 |           | 
#  |    0.7592 |    0.7252 |           | 
#  ------------------|-----------|-----------|-----------|
#  Column Total |    221575 |     18023 |    239598 | 
#  |    0.9248 |    0.0752 |           | 
#  ------------------|-----------|-----------|-----------|
#  Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  104.5371     d.f. =  1     p =  1.542738e-24 
#Pearson's Chi-squared test with Yates' continuity correction 
#------------------------------------------------------------
#Chi^2 =  104.3527     d.f. =  1     p =  1.693268e-24


#total
CrossTable(data$BMI_health,data$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  442399 
#| data$cancer_total_20 
#data$BMI_health |         0 |         1 | Row Total | 
#  ----------------|-----------|-----------|-----------|
#  0 |    100265 |     10242 |    110507 | 
#  |    6.7679 |   72.6704 |           | 
#  |    0.2477 |    0.2717 |           | 
# ----------------|-----------|-----------|-----------|
#  1 |    304443 |     27449 |    331892 | 
# |    2.2534 |   24.1964 |           | 
#  |    0.7523 |    0.7283 |           | 
#  ----------------|-----------|-----------|-----------|
# Column Total |    404708 |     37691 |    442399 | 
# |    0.9148 |    0.0852 |           | 
#  ----------------|-----------|-----------|-----------|
#    Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  105.8881     d.f. =  1     p =  7.801981e-25 
#Pearson's Chi-squared test with Yates' continuity correction 
#------------------------------------------------------------
#Chi^2 =  105.7601     d.f. =  1     p =  8.322475e-25




#############Townsend
##male
CrossTable(male$Townsend_deprivation_index_imp_cat,male$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  202801 
#| male$cancer_total_20 
#male$Townsend_deprivation_index_imp_cat |         0 |         1 | Row Total | 
#  ----------------------------------------|-----------|-----------|-----------|
#  1 |     36677 |      4071 |     40748 | 
#  |    0.3861 |    3.5947 |           | 
#  |    0.2003 |    0.2070 |           | 
#  ----------------------------------------|-----------|-----------|-----------|
#  2 |    108373 |     11700 |    120073 | 
#  |    0.0280 |    0.2608 |           | 
#  |    0.5918 |    0.5949 |           | 
#  ----------------------------------------|-----------|-----------|-----------|
#  3 |     38083 |      3897 |     41980 | 
#  |    0.8014 |    7.4617 |           | 
#  |    0.2080 |    0.1981 |           | 
#  ----------------------------------------|-----------|-----------|-----------|
#  Column Total |    183133 |     19668 |    202801 | 
#  |    0.9030 |    0.0970 |           | 
#  ----------------------------------------|-----------|-----------|-----------|
# Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  12.53257     d.f. =  2     p =  0.001899273 


##female
CrossTable(female$Townsend_deprivation_index_imp_cat,female$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  239598 
#| female$cancer_total_20 
#female$Townsend_deprivation_index_imp_cat |         0 |         1 | Row Total | 
#  ------------------------------------------|-----------|-----------|-----------|
#  1 |     44292 |      3565 |     47857 | 
#  |    0.0275 |    0.3382 |           | 
#  |    0.1999 |    0.1978 |           | 
#  ------------------------------------------|-----------|-----------|-----------|
# 2 |    134249 |     10852 |    145101 | 
#  |    0.0294 |    0.3609 |           | 
#  |    0.6059 |    0.6021 |           | 
#  ------------------------------------------|-----------|-----------|-----------|
#  3 |     43034 |      3606 |     46640 | 
#  |    0.2211 |    2.7182 |           | 
#  |    0.1942 |    0.2001 |           | 
#  ------------------------------------------|-----------|-----------|-----------|
#  Column Total |    221575 |     18023 |    239598 | 
#  |    0.9248 |    0.0752 |           | 
#  ------------------------------------------|-----------|-----------|-----------|
  #  Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  3.695205     d.f. =  2     p =  0.1576146


##Total
CrossTable(data$Townsend_deprivation_index_imp_cat,data$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  442399 
#| data$cancer_total_20 
#data$Townsend_deprivation_index_imp_cat |         0 |         1 | Row Total | 
 # ----------------------------------------|-----------|-----------|-----------|
#  1 |     80969 |      7636 |     88605 | 
#  |    0.0937 |    1.0057 |           | 
#  |    0.2001 |    0.2026 |           | 
#  ----------------------------------------|-----------|-----------|-----------|
#  2 |    242622 |     22552 |    265174 | 
#  |    0.0066 |    0.0708 |           | 
#  |    0.5995 |    0.5983 |           | 
#  ----------------------------------------|-----------|-----------|-----------|
#  3 |     81117 |      7503 |     88620 | 
#  |    0.0274 |    0.2944 |           | 
#  |    0.2004 |    0.1991 |           | 
#  ----------------------------------------|-----------|-----------|-----------|
#  Column Total |    404708 |     37691 |    442399 | 
#  |    0.9148 |    0.0852 |           | 
#  ----------------------------------------|-----------|-----------|-----------|
#   Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#SChi^2 =  1.498579     d.f. =  2     p =  0.4727022



##No current smoking  吸烟状态
#male
CrossTable(male$smoking_health,male$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  202801 
#| male$cancer_total_20 
#male$smoking_health |         0 |         1 | Row Total | 
#  --------------------|-----------|-----------|-----------|
#  0 |     86835 |     10999 |     97834 | 
#  |   25.8391 |  240.5931 |           | 
#  |    0.4742 |    0.5592 |           | 
#  --------------------|-----------|-----------|-----------|
#  1 |     96298 |      8669 |    104967 | 
#  |   24.0832 |  224.2437 |           | 
#  |    0.5258 |    0.4408 |           | 
#  --------------------|-----------|-----------|-----------|
#  Column Total |    183133 |     19668 |    202801 | 
#  |    0.9030 |    0.0970 |           | 
#  --------------------|-----------|-----------|-----------|
  # Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  514.759     d.f. =  1     p =  5.845163e-114 
#Pearson's Chi-squared test with Yates' continuity correction 
#------------------------------------------------------------
#Chi^2 =  514.4184     d.f. =  1     p =  6.932816e-114


#female
CrossTable(female$smoking_health,female$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  239598 
#| female$cancer_total_20 
#female$smoking_health |         0 |         1 | Row Total | 
#  ----------------------|-----------|-----------|-----------|
#  0 |     83660 |      7925 |     91585 | 
#  |   12.6677 |  155.7367 |           | 
#  |    0.3776 |    0.4397 |           | 
#  ----------------------|-----------|-----------|-----------|
#  1 |    137915 |     10098 |    148013 | 
#  |    7.8383 |   96.3641 |           | 
#  |    0.6224 |    0.5603 |           | 
#  ----------------------|-----------|-----------|-----------|
#  Column Total |    221575 |     18023 |    239598 | 
#  |    0.9248 |    0.0752 |           | 
#  ----------------------|-----------|-----------|-----------|
  # Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  272.6068     d.f. =  1     p =  3.06751e-61 
#Pearson's Chi-squared test with Yates' continuity correction 
#------------------------------------------------------------
#Chi^2 =  272.3436     d.f. =  1     p =  3.500497e-61

#total
CrossTable(data$smoking_health,data$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  442399 
#| data$cancer_total_20 
#data$smoking_health |         0 |         1 | Row Total | 
 # --------------------|-----------|-----------|-----------|
  #0 |    170495 |     18924 |    189419 | 
  #|   44.7963 |  481.0008 |           | 
  #|    0.4213 |    0.5021 |           | 
  #--------------------|-----------|-----------|-----------|
  #1 |    234213 |     18767 |    252980 | 
  #|   33.5412 |  360.1498 |           | 
  #|    0.5787 |    0.4979 |           | 
  #--------------------|-----------|-----------|-----------|
  #Column Total |    404708 |     37691 |    442399 | 
  #|    0.9148 |    0.0852 |           | 
#  --------------------|-----------|-----------|-----------|
  # Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  919.4881     d.f. =  1     p =  5.693532e-202 
#Pearson's Chi-squared test with Yates' continuity correction 
#------------------------------------------------------------
#Chi^2 =  919.1581     d.f. =  1     p =  6.716107e-202


####No alcohol consumption  是否饮酒  Alcohol_health
#male
CrossTable(male$Alcohol_health,male$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  202801 
#| male$cancer_total_20 
#male$Alcohol_health |         0 |         1 | Row Total | 
 # --------------------|-----------|-----------|-----------|
#  0 |    177815 |     19267 |    197082 | 
#  |    0.1326 |    1.2350 |           | 
#  |    0.9710 |    0.9796 |           | 
 # --------------------|-----------|-----------|-----------|
#  1 |      5318 |       401 |      5719 | 
#  |    4.5707 |   42.5590 |           | 
#  |    0.0290 |    0.0204 |           | 
#  --------------------|-----------|-----------|-----------|
 # Column Total |    183133 |     19668 |    202801 | 
#  |    0.9030 |    0.0970 |           | 
#  --------------------|-----------|-----------|-----------|
#   Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  48.49734     d.f. =  1     p =  3.307374e-12 
#Pearson's Chi-squared test with Yates' continuity correction 
#------------------------------------------------------------
#Chi^2 =  48.18219     d.f. =  1     p =  3.883993e-12

####No alcohol consumption
#female
CrossTable(female$Alcohol_health,female$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  239598 
#| female$cancer_total_20 
#female$Alcohol_health |         0 |         1 | Row Total | 
 # ----------------------|-----------|-----------|-----------|
#  0 |    208577 |     16989 |    225566 | 
#  |    0.0022 |    0.0273 |           | 
#  |    0.9413 |    0.9426 |           | 
#  ----------------------|-----------|-----------|-----------|
#  1 |     12998 |      1034 |     14032 | 
#  |    0.0357 |    0.4385 |           | 
#  |    0.0587 |    0.0574 |           | 
#  ----------------------|-----------|-----------|-----------|
#  Column Total |    221575 |     18023 |    239598 | 
#  |    0.9248 |    0.0752 |           | 
#  ----------------------|-----------|-----------|-----------|
#   Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  0.5036154     d.f. =  1     p =  0.4779158 
#Pearson's Chi-squared test with Yates' continuity correction 
#------------------------------------------------------------
#Chi^2 =  0.4804774     d.f. =  1     p =  0.4882062


#total
CrossTable(data$Alcohol_health,data$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  442399 
#| data$cancer_total_20 
#data$Alcohol_health |         0 |         1 | Row Total | 
#  --------------------|-----------|-----------|-----------|
#  0 |    386392 |     36256 |    422648 | 
#  |    0.1587 |    1.7042 |           | 
#  |    0.9547 |    0.9619 |           | 
#  --------------------|-----------|-----------|-----------|
#  1 |     18316 |      1435 |     19751 | 
#  |    3.3964 |   36.4687 |           | 
#  |    0.0453 |    0.0381 |           | 
#  --------------------|-----------|-----------|-----------|
#  Column Total |    404708 |     37691 |    442399 | 
#  |    0.9148 |    0.0852 |           | 
 # --------------------|-----------|-----------|-----------|
#  Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  41.72799     d.f. =  1     p =  1.048965e-10 
#Pearson's Chi-squared test with Yates' continuity correction 
#------------------------------------------------------------
#Chi^2 =  41.55971     d.f. =  1     p =  1.14325e-10




####Regular physical activity
#male
CrossTable(male$physical_health,male$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  202801 
#| male$cancer_total_20 
#male$physical_health |         0 |         1 | Row Total | 
#  ---------------------|-----------|-----------|-----------|
#  0 |     64644 |      7141 |     71785 | 
#  |    0.4952 |    4.6108 |           | 
#  |    0.3530 |    0.3631 |           | 
#  ---------------------|-----------|-----------|-----------|
#  1 |    118489 |     12527 |    131016 | 
#  |    0.2713 |    2.5263 |           | 
#  |    0.6470 |    0.6369 |           | 
#  ---------------------|-----------|-----------|-----------|
#  Column Total |    183133 |     19668 |    202801 | 
#  |    0.9030 |    0.0970 |           | 
#  ---------------------|-----------|-----------|-----------|
#  Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  7.903602     d.f. =  1     p =  0.004933645 
#Pearson's Chi-squared test with Yates' continuity correction 
#------------------------------------------------------------
#Chi^2 =  7.85955     d.f. =  1     p =  0.005055297


#female
CrossTable(female$physical_health,female$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  239598 
#| female$cancer_total_20 
#female$physical_health |         0 |         1 | Row Total | 
#  -----------------------|-----------|-----------|-----------|
#  0 |     81266 |      6826 |     88092 | 
#  |    0.4888 |    6.0098 |           | 
#  |    0.3668 |    0.3787 |           | 
#  -----------------------|-----------|-----------|-----------|
#  1 |    140309 |     11197 |    151506 | 
#  |    0.2842 |    3.4944 |           | 
#  |    0.6332 |    0.6213 |           | 
#  -----------------------|-----------|-----------|-----------|
#  Column Total |    221575 |     18023 |    239598 | 
#  |    0.9248 |    0.0752 |           | 
#  -----------------------|-----------|-----------|-----------|
#  Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  10.27723     d.f. =  1     p =  0.001346819 
#Pearson's Chi-squared test with Yates' continuity correction 
#------------------------------------------------------------
#Chi^2 =  10.2258     d.f. =  1     p =  0.001384899 


#total
CrossTable(data$physical_health,data$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  442399 
#| data$cancer_total_20 
#data$physical_health |         0 |         1 | Row Total | 
#  ---------------------|-----------|-----------|-----------|
#  0 |    145910 |     13967 |    159877 | 
#  |    0.8185 |    8.7882 |           | 
#  |    0.3605 |    0.3706 |           | 
#  ---------------------|-----------|-----------|-----------|
#  1 |    258798 |     23724 |    282522 | 
#  |    0.4632 |    4.9732 |           | 
#  |    0.6395 |    0.6294 |           | 
#  ---------------------|-----------|-----------|-----------|
#  Column Total |    404708 |     37691 |    442399 | 
#  |    0.9148 |    0.0852 |           | 
#  ---------------------|-----------|-----------|-----------|
#  Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  15.04305     d.f. =  1     p =  0.0001050865 
#Pearson's Chi-squared test with Yates' continuity correction 
#------------------------------------------------------------
#Chi^2 =  14.9996     d.f. =  1     p =  0.0001075339



####家族史数量等级分布####

table(male_inci$FH_total_2)
#0     1 
#13255  6413 

#6413/19668   0.3261


table(male_con$FH_total_2)
#0      1 
#133269  49864 
#49864/183133     0.2723


table(female_inci$FH_total_2)
#0     1 
#11737  6286 
#6286/18023   0.3488

table(female_con$FH_total_2)
#0      1 
#155797  65778
#65778/221575    0.2969


table(inci$FH_total_2)
#0     1 
#24992 12699 

#12699/37691    0.3369  


table(con$FH_total_2)
#0      1 
#289066 115642
#115642/404708   0.2857



#男性  肠癌  肺癌  前列腺癌  三分类 0，1，≥2
#肠癌
table(male_inci$FH_Bowel_3)
#0      1      2 
#17100  2398   170 
table(male_con$FH_Bowel_3)
#162948  18988   1197


#肺癌
table(male_inci$FH_Lung_3)
#0      1      2  
#16915  2542   211
table(male_con$FH_Lung_3)
#161231  20370   1532

#前列腺癌
table(male_inci$FH_Prostate_3)
#0      1      2 
#17672  1910    86
table(male_con$FH_Prostate_3)
#169356  13468    309

#女性  肠癌  肺癌  乳腺癌  三分类 0，1，≥2
#肠癌
table(female_inci$FH_Bowel_3)
#0      1      2 
#15760  2134   129 
table(female_con$FH_Bowel_3)
#198040  22292   1243

#肺癌
table(female_inci$FH_Lung_3)
#0      1      2 
#15426  2393   204 
table(female_con$FH_Lung_3)
#193823  26009   1743

#乳腺癌
table(female_inci$FH_Breast_3)
#0      1      2 
#15562  2319   142 
table(female_con$FH_Breast_3)
#197811  22744   1020


#总体  肠癌  肺癌   三分类 0，1，≥2
#肠癌
table(inci$FH_Bowel_3)
#0      1      2 
#32860  4532   299 
table(con$FH_Bowel_3)
#360988  41280   2440


#肺癌
table(inci$FH_Lung_3)
#0      1      2 
#32341  4935   415 
table(con$FH_Lung_3)
#355054  46379   3275



####多肿瘤家族史  四分类 0，1，2，≥3####
#男性
table(male_inci$FH_total_4)
#0      1      2      3 
#13255  5220  1033   160 
table(male_con$FH_total_4)
#133269  41867   6993   1004

#女性
table(female_inci$FH_total_4)
#0      1      2      3 
#11737  4975  1114   197 
table(female_con$FH_total_4)
#155797  54113  10167   1498

#总体
table(inci$FH_total_4)
#0      1      2      3 
#24992 10195  2147   357
table(con$FH_total_4)
#289066  95980  17160   2502






####PRS   #colorectal: COL_PRS lung: LC_PRS  prostate:  PRO_PRS  breast:  BC_PRS####

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3

#男性（case/control）-肠癌 肺癌 前列腺癌
#肠癌
Q20=quantile(male$COL_PRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(male$COL_PRS,seq(0.05,1,0.05))[[16]]

male$COL_PRS_risk_3_Q20=0
male[male$COL_PRS<=Q20,]$COL_PRS_risk_3_Q20=1
male[male$COL_PRS>Q20&male$COL_PRS<Q80,]$COL_PRS_risk_3_Q20=2
male[male$COL_PRS>=Q80,]$COL_PRS_risk_3_Q20=3
male$COL_PRS_risk_3_Q20f=as.factor(male$COL_PRS_risk_3_Q20)
table(male$COL_PRS_risk_3_Q20f)
#1      2      3 
#40561 121679  40561 

round(mean(male_inci$COL_PRS),2)#3.22
round(sd(male_inci$COL_PRS),2)#0.4

round(mean(male_con$COL_PRS),2)#3.18
round(sd(male_con$COL_PRS),2)#0.39

CrossTable(male$COL_PRS_risk_3_Q20,male$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)

#Total Observations in Table:  202801 

#                        | male$cancer_total_20 
#male$COL_PRS_risk_3_Q20 |         0 |         1 | Row Total | 
#------------------------|-----------|-----------|-----------|
#1                       |     37036 |      3525 |     40561 | 
#                        |    4.5599 |   42.4583 |           | 
#                        |    0.2022 |    0.1792 |           | 
#------------------------|-----------|-----------|-----------|
#2                       |    109974 |     11705 |    121679 | 
#                        |    0.0833 |    0.7752 |           | 
#                        |    0.6005 |    0.5951 |           | 
#------------------------|-----------|-----------|-----------|
#3                       |     36123 |      4438 |     40561 | 
#                        |    6.9440 |   64.6573 |           | 
#                        |    0.1973 |    0.2256 |           | 
#------------------------|-----------|-----------|-----------|
#           Column Total |    183133 |     19668 |    202801 | 
#                        |    0.9030 |    0.0970 |           | 
#------------------------|-----------|-----------|-----------|
  
#    Statistics for All Table Factors
#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  119.4781     d.f. =  2     p =  1.136761e-26


#肺癌  LC_PRS
Q20=quantile(male$LC_PRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(male$LC_PRS,seq(0.05,1,0.05))[[16]]

male$LC_PRS_risk_3_Q20=0
male[male$LC_PRS<=Q20,]$LC_PRS_risk_3_Q20=1
male[male$LC_PRS>Q20&male$LC_PRS<Q80,]$LC_PRS_risk_3_Q20=2
male[male$LC_PRS>=Q80,]$LC_PRS_risk_3_Q20=3
male$LC_PRS_risk_3_Q20f=as.factor(male$LC_PRS_risk_3_Q20)
table(male$LC_PRS_risk_3_Q20f)
#1      2      3 
#40561 121679  40561  

round(mean(male_inci$LC_PRS),2)#0.09
round(sd(male_inci$LC_PRS),2)#0.11

round(mean(male_con$LC_PRS),2)#0.08
round(sd(male_con$LC_PRS),2)# 0.11

CrossTable(male$LC_PRS_risk_3_Q20,male$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  202801 

#                       | male$cancer_total_20 
#male$LC_PRS_risk_3_Q20 |         0 |         1 | Row Total | 
#-----------------------|-----------|-----------|-----------|
#                     1 |     36730 |      3831 |     40561 | 
#                       |    0.2878 |    2.6801 |           | 
#                       |    0.2006 |    0.1948 |           | 
#-----------------------|-----------|-----------|-----------|
#                     2 |    109875 |     11804 |    121679 | 
#                       |    0.0001 |    0.0010 |           | 
#                       |    0.6000 |    0.6002 |           | 
#-----------------------|-----------|-----------|-----------|
#                     3 |     36528 |      4033 |     40561 | 
#                       |    0.2693 |    2.5078 |           | 
#                       |    0.1995 |    0.2051 |           | 
#-----------------------|-----------|-----------|-----------|
#          Column Total |    183133 |     19668 |    202801 | 
#                       |    0.9030 |    0.0970 |           | 
#-----------------------|-----------|-----------|-----------|
  
#  Statistics for All Table Factors

#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  5.746152     d.f. =  2     p =  0.05652479 


#前列腺癌  PRO_PRS
Q20=quantile(male$PRO_PRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(male$PRO_PRS,seq(0.05,1,0.05))[[16]]

male$PRO_PRS_risk_3_Q20=0
male[male$PRO_PRS<=Q20,]$PRO_PRS_risk_3_Q20=1
male[male$PRO_PRS>Q20&male$PRO_PRS<Q80,]$PRO_PRS_risk_3_Q20=2
male[male$PRO_PRS>=Q80,]$PRO_PRS_risk_3_Q20=3
male$PRO_PRS_risk_3_Q20f=as.factor(male$PRO_PRS_risk_3_Q20)
table(male$PRO_PRS_risk_3_Q20f)
#1      2      3 
#40561 121679  40561  

round(mean(male_inci$PRO_PRS),2)#8.41
round(sd(male_con$PRO_PRS),2)#0.69

round(mean(male_con$PRO_PRS),2)#8.20
round(sd(male_con$PRO_PRS),2)#0.69

CrossTable(male$PRO_PRS_risk_3_Q20,male$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)

#Total Observations in Table:  202801 
#                        | male$cancer_total_20 
#male$PRO_PRS_risk_3_Q20 |         0 |         1 | Row Total | 
#------------------------|-----------|-----------|-----------|
#                      1 |     37740 |      2821 |     40561 | 
#                        |   33.8013 |  314.7313 |           | 
#                        |    0.2061 |    0.1434 |           | 
#------------------------|-----------|-----------|-----------|
#                      2 |    110503 |     11176 |    121679 | 
#                        |    3.5510 |   33.0644 |           | 
#                        |    0.6034 |    0.5682 |           | 
#------------------------|-----------|-----------|-----------|
#                      3 |     34890 |      5671 |     40561 | 
#                        |   82.4054 |  767.2945 |           | 
#                        |    0.1905 |    0.2883 |           | 
#------------------------|-----------|-----------|-----------|
#          Column Total  |    183133 |     19668 |    202801 | 
#                        |    0.9030 |    0.0970 |           | 
#------------------------|-----------|-----------|-----------|
  
  
#  Statistics for All Table Factors


#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  1234.848     d.f. =  2     p =  7.180951e-269 



#CPRS
round(mean(male_inci$CPRS),2)#3.80
round(sd(male_inci$CPRS),2)#0.21

round(mean(male_con$CPRS),2)#3.74
round(sd(male_con$CPRS),2)#0.20

Q20=quantile(male$CPRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(male$CPRS,seq(0.05,1,0.05))[[16]]
male$CPRS_risk_3_Q20=0
male[male$CPRS<=Q20,]$CPRS_risk_3_Q20=1
male[male$CPRS>Q20&male$CPRS<Q80,]$CPRS_risk_3_Q20=2
male[male$CPRS>=Q80,]$CPRS_risk_3_Q20=3
male$CPRS_risk_3_Q20f=as.factor(male$CPRS_risk_3_Q20)
table(male$CPRS_risk_3_Q20)
#1      2      3 
#40561 121679  40561 

CrossTable(male$CPRS_risk_3_Q20,male$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  202801 

#                     | male$cancer_total_20 
#male$CPRS_risk_3_Q20 |         0 |         1 | Row Total | 
#---------------------|-----------|-----------|-----------|
#                   1 |     37914 |      2647 |     40561 | 
#                     |   45.1996 |  420.8630 |           | 
#                     |    0.2070 |    0.1346 |           | 
#---------------------|-----------|-----------|-----------|
#                   2 |    110463 |     11216 |    121679 | 
#                     |    3.1108 |   28.9653 |           | 
#                     |    0.6032 |    0.5703 |           | 
#---------------------|-----------|-----------|-----------|
#                   3 |     34756 |      5805 |     40561 | 
#                     |   95.6075 |  890.2223 |           | 
#                     |    0.1898 |    0.2951 |           | 
#---------------------|-----------|-----------|-----------|
#        Column Total |    183133 |     19668 |    202801 | 
#                     |    0.9030 |    0.0970 |           | 
#---------------------|-----------|-----------|-----------|
  
  
#  Statistics for All Table Factors


#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  1483.969     d.f. =  2     p =  5.928788e-323



####女性（case/control）-肠癌 肺癌 乳腺癌####
#肠癌
Q20=quantile(female$COL_PRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(female$COL_PRS,seq(0.05,1,0.05))[[16]]

female$COL_PRS_risk_3_Q20=0
female[female$COL_PRS<=Q20,]$COL_PRS_risk_3_Q20=1
female[female$COL_PRS>Q20&female$COL_PRS<Q80,]$COL_PRS_risk_3_Q20=2
female[female$COL_PRS>=Q80,]$COL_PRS_risk_3_Q20=3
female$COL_PRS_risk_3_Q20f=as.factor(female$COL_PRS_risk_3_Q20)
table(female$COL_PRS_risk_3_Q20f)
#1      2      3 
#47920 143758  47920

round(mean(female_inci$COL_PRS),2)#3.21
round(sd(female_inci$COL_PRS),2)#0.40

round(mean(female_con$COL_PRS),2)#3.19
round(sd(female_con$COL_PRS),2)#0.39

CrossTable(female$COL_PRS_risk_3_Q20,female$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)

#Total Observations in Table:  239598 
#                          | female$cancer_total_20 
#female$COL_PRS_risk_3_Q20 |         0 |         1 | Row Total | 
#--------------------------|-----------|-----------|-----------|
#                        1 |     44589 |      3331 |     47920 | 
#                          |    1.6896 |   20.7715 |           | 
#                          |    0.2012 |    0.1848 |           | 
#--------------------------|-----------|-----------|-----------|
#                        2 |    133027 |     10731 |    143758 | 
#                          |    0.0515 |    0.6331 |           | 
#                          |    0.6004 |    0.5954 |           | 
#--------------------------|-----------|-----------|-----------|
#                        3 |     43959 |      3961 |     47920 | 
#                          |    2.8658 |   35.2323 |           | 
#                          |    0.1984 |    0.2198 |           | 
#--------------------------|-----------|-----------|-----------|
#             Column Total |    221575 |     18023 |    239598 | 
#                          |    0.9248 |    0.0752 |           | 
#--------------------------|-----------|-----------|-----------|
  
  
#  Statistics for All Table Factors


#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  61.24373     d.f. =  2     p =  5.024514e-14 


#肺癌  LC_PRS
Q20=quantile(female$LC_PRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(female$LC_PRS,seq(0.05,1,0.05))[[16]]

female$LC_PRS_risk_3_Q20=0
female[female$LC_PRS<=Q20,]$LC_PRS_risk_3_Q20=1
female[female$LC_PRS>Q20&female$LC_PRS<Q80,]$LC_PRS_risk_3_Q20=2
female[female$LC_PRS>=Q80,]$LC_PRS_risk_3_Q20=3
female$LC_PRS_risk_3_Q20f=as.factor(female$LC_PRS_risk_3_Q20)
table(female$LC_PRS_risk_3_Q20f)
#1      2      3 
#47920 143758  47920  

round(mean(female_inci$LC_PRS),2)#0.09
round(sd(female_inci$LC_PRS),2)#0.11

round(mean(female_con$LC_PRS),2)#0.09
round(sd(female_con$LC_PRS),2)#0.11

CrossTable(female$LC_PRS_risk_3_Q20,female$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  239598 
#                         | female$cancer_total_20 
#female$LC_PRS_risk_3_Q20 |         0 |         1 | Row Total | 
#-------------------------|-----------|-----------|-----------|
#                       1 |     44392 |      3528 |     47920 | 
#                         |    0.1325 |    1.6291 |           | 
#                         |    0.2003 |    0.1957 |           | 
#-------------------------|-----------|-----------|-----------|
#                       2 |    133095 |     10663 |    143758 | 
#                         |    0.1709 |    2.1013 |           | 
#                         |    0.6007 |    0.5916 |           | 
#-------------------------|-----------|-----------|-----------|
#                       3 |     44088 |      3832 |     47920 | 
#                         |    1.1666 |   14.3419 |           | 
#                         |    0.1990 |    0.2126 |           | 
#-------------------------|-----------|-----------|-----------|
#            Column Total |    221575 |     18023 |    239598 | 
#                         |    0.9248 |    0.0752 |           | 
#-------------------------|-----------|-----------|-----------|
  
  
#  Statistics for All Table Factors


#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  19.54218     d.f. =  2     p =  5.707821e-05 


#乳腺癌  BC_PRS
Q20=quantile(female$BC_PRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(female$BC_PRS,seq(0.05,1,0.05))[[16]]

female$BC_PRS_risk_3_Q20=0
female[female$BC_PRS<=Q20,]$BC_PRS_risk_3_Q20=1
female[female$BC_PRS>Q20&female$BC_PRS<Q80,]$BC_PRS_risk_3_Q20=2
female[female$BC_PRS>=Q80,]$BC_PRS_risk_3_Q20=3
female$BC_PRS_risk_3_Q20f=as.factor(female$BC_PRS_risk_3_Q20)
table(female$BC_PRS_risk_3_Q20f)
#1      2      3 
#47920 143758  47920  

round(mean(female_inci$BC_PRS),2)#-0.09
round(sd(female_inci$BC_PRS),2)#0.13

round(mean(female_con$BC_PRS),2)#-0.1
round(sd(female_con$BC_PRS),2)#0.13

CrossTable(female$BC_PRS_risk_3_Q20,female$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  239598 
#                         | female$cancer_total_20 
#female$BC_PRS_risk_3_Q20 |         0 |         1 | Row Total | 
#-------------------------|-----------|-----------|-----------|
#                       1 |     44995 |      2925 |     47920 | 
#                         |   10.4230 |  128.1399 |           | 
#                         |    0.2031 |    0.1623 |           | 
#-------------------------|-----------|-----------|-----------|
#                       2 |    132983 |     10775 |    143758 | 
#                         |    0.0113 |    0.1388 |           | 
#                         |    0.6002 |    0.5978 |           | 
#-------------------------|-----------|-----------|-----------|
#                       3 |     43597 |      4323 |     47920 | 
#                         |   11.6451 |  143.1646 |           | 
#                         |    0.1968 |    0.2399 |           | 
#-------------------------|-----------|-----------|-----------|
#            Column Total |    221575 |     18023 |    239598 | 
#                         |    0.9248 |    0.0752 |           | 
#-------------------------|-----------|-----------|-----------|
  
  
#  Statistics for All Table Factors


#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  293.5226     d.f. =  2     p =  1.829691e-64 


#CPRS
round(mean(female_inci$CPRS),2)#1.45
round(sd(female_inci$CPRS),2)#0.08

round(mean(female_con$CPRS),2)#1.44
round(sd(female_con$CPRS),2)#0.08

Q20=quantile(female$CPRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(female$CPRS,seq(0.05,1,0.05))[[16]]

female$CPRS_risk_3_Q20=0
female[female$CPRS<=Q20,]$CPRS_risk_3_Q20=1
female[female$CPRS>Q20&female$CPRS<Q80,]$CPRS_risk_3_Q20=2
female[female$CPRS>=Q80,]$CPRS_risk_3_Q20=3
female$CPRS_risk_3_Q20f=as.factor(female$CPRS_risk_3_Q20)
table(female$CPRS_risk_3_Q20)
#1      2      3 
#47920 143758  47920  

CrossTable(female$CPRS_risk_3_Q20,female$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  239598 
#                       | female$cancer_total_20 
#female$CPRS_risk_3_Q20 |         0 |         1 | Row Total | 
#-----------------------|-----------|-----------|-----------|
#                     1 |     45095 |      2825 |     47920 | 
#                       |   13.7159 |  168.6229 |           | 
#                       |    0.2035 |    0.1567 |           | 
#-----------------------|-----------|-----------|-----------|
#                     2 |    133106 |     10652 |    143758 | 
#                       |    0.1968 |    2.4191 |           | 
#                       |    0.6007 |    0.5910 |           | 
#-----------------------|-----------|-----------|-----------|
#                     3 |     43374 |      4546 |     47920 | 
#                       |   19.9971 |  245.8442 |           | 
#                       |    0.1958 |    0.2522 |           | 
#-----------------------|-----------|-----------|-----------|
#          Column Total |    221575 |     18023 |    239598 | 
#                       |    0.9248 |    0.0752 |           | 
#-----------------------|-----------|-----------|-----------|
  
  
#  Statistics for All Table Factors


#Pearson's Chi-squared test 
------------------------------------------------------------
#Chi^2 =  450.7959     d.f. =  2     p =  1.29099e-98 



####总体（case/control）-肠癌 肺癌  CPRS####
#肠癌
Q20=quantile(data$COL_PRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(data$COL_PRS,seq(0.05,1,0.05))[[16]]

data$COL_PRS_risk_3_Q20=0
data[data$COL_PRS<=Q20,]$COL_PRS_risk_3_Q20=1
data[data$COL_PRS>Q20&data$COL_PRS<Q80,]$COL_PRS_risk_3_Q20=2
data[data$COL_PRS>=Q80,]$COL_PRS_risk_3_Q20=3
data$COL_PRS_risk_3_Q20f=as.factor(data$COL_PRS_risk_3_Q20)
table(data$COL_PRS_risk_3_Q20f)
#1      2      3 
#88480 265439  88480

round(mean(inci$COL_PRS),2)#3.22
round(sd(inci$COL_PRS),2)#0.4

round(mean(con$COL_PRS),2)#3.18
round(sd(con$COL_PRS),2)#0.39

CrossTable(data$COL_PRS_risk_3_Q20,data$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  442399 
##                       | data$cancer_total_20 
#data$COL_PRS_risk_3_Q20 |         0 |         1 | Row Total | 
#------------------------|-----------|-----------|-----------|
#                      1 |     81624 |      6856 |     88480 | 
#                        |    5.7501 |   61.7414 |           | 
#                        |    0.2017 |    0.1819 |           | 
#------------------------|-----------|-----------|-----------|
#                      2 |    243005 |     22434 |    265439 | 
#                        |    0.1343 |    1.4417 |           | 
#                        |    0.6004 |    0.5952 |           | 
#------------------------|-----------|-----------|-----------|
#                      3 |     80079 |      8401 |     88480 | 
#                        |    9.1967 |   98.7494 |           | 
#                        |    0.1979 |    0.2229 |           | 
#------------------------|-----------|-----------|-----------|
#           Column Total |    404708 |     37691 |    442399 | 
#                        |    0.9148 |    0.0852 |           | 
#------------------------|-----------|-----------|-----------|
  
  
#  Statistics for All Table Factors


#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  177.0135     d.f. =  2     p =  3.647546e-39 
  
  
  
#肺癌  LC_PRS
Q20=quantile(data$LC_PRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(data$LC_PRS,seq(0.05,1,0.05))[[16]]

data$LC_PRS_risk_3_Q20=0
data[data$LC_PRS<=Q20,]$LC_PRS_risk_3_Q20=1
data[data$LC_PRS>Q20&data$LC_PRS<Q80,]$LC_PRS_risk_3_Q20=2
data[data$LC_PRS>=Q80,]$LC_PRS_risk_3_Q20=3
data$LC_PRS_risk_3_Q20f=as.factor(data$LC_PRS_risk_3_Q20)
table(data$LC_PRS_risk_3_Q20f)
#1      2      3 
#88480 265439  88480 

round(mean(inci$LC_PRS),2)#0.09
round(sd(inci$LC_PRS),2)#0.11

round(mean(con$LC_PRS),2)#0.08
round(sd(con$LC_PRS),2)#0.11

CrossTable(data$LC_PRS_risk_3_Q20,data$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  442399 
#                       | data$cancer_total_20 
#data$LC_PRS_risk_3_Q20 |         0 |         1 | Row Total | 
#-----------------------|-----------|-----------|-----------|
#                     1 |     81121 |      7359 |     88480 | 
#                       |    0.3968 |    4.2608 |           | 
#                       |    0.2004 |    0.1952 |           | 
#-----------------------|-----------|-----------|-----------|
#                     2 |    242977 |     22462 |    265439 | 
#                       |    0.0959 |    1.0293 |           | 
#                       |    0.6004 |    0.5960 |           | 
#-----------------------|-----------|-----------|-----------|
#                     3 |     80610 |      7870 |     88480 | 
#                       |    1.3600 |   14.6029 |           | 
#                       |    0.1992 |    0.2088 |           | 
#-----------------------|-----------|-----------|-----------|
#          Column Total |    404708 |     37691 |    442399 | 
#                       |    0.9148 |    0.0852 |           | 
#-----------------------|-----------|-----------|-----------|
  
  
#  Statistics for All Table Factors


#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  21.74562     d.f. =  2     p =  1.896695e-05  



#CPRS
round(mean(inci$CPRS),2)#2.68
round(sd(inci$CPRS),2)# 1.18

round(mean(con$CPRS),2)#2.48
round(sd(con$CPRS),2)#1.15

Q20=quantile(data$CPRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(data$CPRS,seq(0.05,1,0.05))[[16]]

data$CPRS_risk_3_Q20=0
data[data$CPRS<=Q20,]$CPRS_risk_3_Q20=1
data[data$CPRS>Q20&data$CPRS<Q80,]$CPRS_risk_3_Q20=2
data[data$CPRS>=Q80,]$CPRS_risk_3_Q20=3
data$CPRS_risk_3_Q20f=as.factor(data$CPRS_risk_3_Q20)
table(data$CPRS_risk_3_Q20)
#1      2      3 
#88480 265439  88480 

CrossTable(data$CPRS_risk_3_Q20,data$cancer_total_20,digit=4,prop.r=F,prop.c=T,prop.t=F,chisq=T)
#Total Observations in Table:  442399 
#                     | data$cancer_total_20 
#data$CPRS_risk_3_Q20 |         0 |         1 | Row Total | 
#---------------------|-----------|-----------|-----------|
#                   1 |     82954 |      5526 |     88480 | 
#                     |   50.0238 |  537.1320 |           | 
#                     |    0.2050 |    0.1466 |           | 
#---------------------|-----------|-----------|-----------|
#                   2 |    244178 |     21261 |    265439 | 
#                     |    7.5451 |   81.0160 |           | 
#                     |    0.6033 |    0.5641 |           | 
#---------------------|-----------|-----------|-----------|
#                   3 |     77576 |     10904 |     88480 | 
#                     |  139.9586 | 1502.8083 |           | 
#                     |    0.1917 |    0.2893 |           | 
#---------------------|-----------|-----------|-----------|
#        Column Total |    404708 |     37691 |    442399 | 
#         |    0.9148 |    0.0852 |           | 
#---------------------|-----------|-----------|-----------|
  
  
#  Statistics for All Table Factors


#Pearson's Chi-squared test 
#------------------------------------------------------------
#Chi^2 =  2318.484     d.f. =  2     p =  0
