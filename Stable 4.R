#Stable 4
#PRS与对应肿瘤发生风险


#加载数据
load("imputed_final_20240120_UKB_database_for_analysis_202801_without_prevalent_cancer_male.rdata")
load("imputed_final_20240120_UKB_database_for_analysis_239598_without_prevalent_cancer_female.rdata")


#加载R包
library("survival")
library("survminer")
library("rstpm2")
library("flexsurv")
library("rstpm2")
library("rms")


#肠癌 COL_PRS #

#male#
Q20=quantile(male$COL_PRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(male$COL_PRS,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
male$COL_PRS_3_Q20=0
male[male$COL_PRS<=Q20,]$COL_PRS_3_Q20=1
male[male$COL_PRS>Q20&male$COL_PRS<Q80,]$COL_PRS_3_Q20=2
male[male$COL_PRS>=Q80,]$COL_PRS_3_Q20=3
table(male$COL_PRS_3_Q20)
#1      2      3 
#40561 121679  40561 

#人年
fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male[male$COL_PRS_3_Q20==1,],scale = 1)
fit2$pyears#
fit2$event#
#285/437937

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male[male$COL_PRS_3_Q20==2,],scale = 1)
fit2$pyears#
fit2$event#
#1634/1311358

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=male[male$COL_PRS_3_Q20==3,],scale = 1)
fit2$pyears#
fit2$event#
#919/434792

#HR
model11 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(COL_PRS_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model11)
res <- summary(model11)
print(res$coefficients,digits = 2)
print(res$conf.int,digit=2)
#                                        coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(COL_PRS_3_Q20)2            0.63060      1.88  0.06427  9.811  1.0e-22
#as.factor(COL_PRS_3_Q20)3            1.14957      3.16  0.06800 16.905  4.2e-64

#exp(coef) exp(-coef) lower .95 upper .95
#as.factor(COL_PRS_3_Q20)2                1.88       0.53      1.66       2.13
#as.factor(COL_PRS_3_Q20)3                3.16       0.32      2.76       3.61

#p for trend
model12 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.numeric(COL_PRS_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
res <- summary(model12)
print(res$coefficients,digits = 2)
print(res$conf.int,digit=2)
#                                        coef exp(coef) se(coef)      z Pr(>|z|)
#as.numeric(COL_PRS_3_Q20)            0.55603      1.74  0.03033 18.335  4.4e-75
#                                    exp(coef) exp(-coef) lower .95 upper .95
#as.numeric(COL_PRS_3_Q20)              1.74       0.57      1.64       1.9



#female#
Q20=quantile(female$COL_PRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(female$COL_PRS,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
female$COL_PRS_3_Q20=0
female[female$COL_PRS<=Q20,]$COL_PRS_3_Q20=1
female[female$COL_PRS>Q20&female$COL_PRS<Q80,]$COL_PRS_3_Q20=2
female[female$COL_PRS>=Q80,]$COL_PRS_3_Q20=3
table(female$COL_PRS_3_Q20)
#1      2      3 
#47920 143758  47920 

#人年
fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female[female$COL_PRS_3_Q20==1,],scale = 1)
fit2$pyears#
fit2$event#
#221/524152

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female[female$COL_PRS_3_Q20==2,],scale = 1)
fit2$pyears#
fit2$event#
#1219/1573996

fit2 <- pyears(Surv(CRC_difftime_new,CRC_total.y)~1,data=female[female$COL_PRS_3_Q20==3,],scale = 1)
fit2$pyears#
fit2$event#
#727/523542

#HR
model11 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.factor(COL_PRS_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
res <- summary(model11)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                       coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(COL_PRS_3_Q20)2            6.0e-01      1.82  0.07324  8.209  2.2e-16
#as.factor(COL_PRS_3_Q20)3            1.2e+00      3.27  0.07716 15.350  3.5e-53

#                                       exp(coef) exp(-coef) lower .95 upper .95
#as.factor(COL_PRS_3_Q20)2                1.82       0.55      1.58       2.11
#as.factor(COL_PRS_3_Q20)3                3.27       0.31      2.81       3.80

#p for trend
model12 <- coxph(Surv(CRC_difftime_new,CRC_total.y)~as.numeric(COL_PRS_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
res=summary(model12)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                        coef exp(coef) se(coef)      z Pr(>|z|)
#as.numeric(COL_PRS_3_Q20)            5.9e-01      1.80  0.03480 16.928  2.8e-64
#                                        exp(coef) exp(-coef) lower .95 upper .95
#as.numeric(COL_PRS_3_Q20)                1.80       0.55      1.68       1.9





#肺癌   LC_PRS#

#male#
Q20=quantile(male$LC_PRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(male$LC_PRS,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
male$LC_PRS_3_Q20=0
male[male$LC_PRS<=Q20,]$LC_PRS_3_Q20=1
male[male$LC_PRS>Q20&male$LC_PRS<Q80,]$LC_PRS_3_Q20=2
male[male$LC_PRS>=Q80,]$LC_PRS_3_Q20=3
table(male$LC_PRS_3_Q20)
#1      2      3 
#40576 121664  40561 

#人年
fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male[male$LC_PRS_3_Q20==1,],scale = 1)
fit2$pyears#
fit2$event#
#226/438915

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male[male$LC_PRS_3_Q20==2,],scale = 1)
fit2$pyears#
fit2$event#
#1097/1315520

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=male[male$LC_PRS_3_Q20==3,],scale = 1)
fit2$pyears#
fit2$event#
#484/438326


#HR
model11 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(LC_PRS_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
summary(model11)
res=summary(model11)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                        coef exp(coef) se(coef)     z Pr(>|z|)
#as.factor(LC_PRS_3_Q20)2             0.50391      1.66  0.07306  6.90  5.3e-12
#as.factor(LC_PRS_3_Q20)3             0.78912      2.20  0.08064  9.79  1.3e-22

#exp(coef) exp(-coef) lower .95 upper .95
#as.factor(LC_PRS_3_Q20)2                 1.66       0.60      1.43       1.91
#as.factor(LC_PRS_3_Q20)3                 2.20       0.45      1.88       2.58

#p for trend
model12 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.numeric(LC_PRS_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
res=summary(model12)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#coef exp(coef) se(coef)     z Pr(>|z|)
#as.numeric(LC_PRS_3_Q20)             0.37048      1.45  0.03740  9.91  3.9e-23

#exp(coef) exp(-coef) lower .95 upper .95
#as.numeric(LC_PRS_3_Q20)                 1.45       0.69      1.35       1.6




#female#
Q20=quantile(female$LC_PRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(female$LC_PRS,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
female$LC_PRS_3_Q20=0
female[female$LC_PRS<=Q20,]$LC_PRS_3_Q20=1
female[female$LC_PRS>Q20&female$LC_PRS<Q80,]$LC_PRS_3_Q20=2
female[female$LC_PRS>=Q80,]$LC_PRS_3_Q20=3
table(female$LC_PRS_3_Q20)
#1      2      3 
#47920 143758  47920 

#人年
fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female[female$LC_PRS_3_Q20==1,],scale = 1)
fit2$pyears#
fit2$event#
#231/525960


fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female[female$LC_PRS_3_Q20==2,],scale = 1)
fit2$pyears#
fit2$event#
#1020/1576184

fit2 <- pyears(Surv(LC_difftime_new,LC_total.y)~1,data=female[female$LC_PRS_3_Q20==3,],scale = 1)
fit2$pyears#
fit2$event#
#431/525235


#HR
model11 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.factor(LC_PRS_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
res=summary(model11)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                        coef exp(coef) se(coef)     z Pr(>|z|)
#as.factor(LC_PRS_3_Q20)2             0.38958      1.48  0.07288  5.35  9.0e-08
#as.factor(LC_PRS_3_Q20)3             0.63429      1.89  0.08159  7.77  7.6e-15

#                                        exp(coef) exp(-coef) lower .95 upper .95
#as.factor(LC_PRS_3_Q20)2                 1.48       0.68      1.28       1.70
#as.factor(LC_PRS_3_Q20)3                 1.89       0.53      1.61       2.21

#p for trend
model12 <- coxph(Surv(LC_difftime_new,LC_total.y)~as.numeric(LC_PRS_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
res=summary(model12)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                     coef exp(coef) se(coef)      z Pr(>|z|)
#as.numeric(LC_PRS_3_Q20)             0.3039      1.36  0.03878  7.835  4.7e-15

#                                     exp(coef) exp(-coef) lower .95 upper .95
#as.numeric(LC_PRS_3_Q20)                 1.36       0.74      1.26       1.5





#前列腺癌 PRO_PRS#

#male#
Q20=quantile(male$PRO_PRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(male$PRO_PRS,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
male$PRO_PRS_3_Q20=0
male[male$PRO_PRS<=Q20,]$PRO_PRS_3_Q20=1
male[male$PRO_PRS>Q20&male$PRO_PRS<Q80,]$PRO_PRS_3_Q20=2
male[male$PRO_PRS>=Q80,]$PRO_PRS_3_Q20=3
table(male$PRO_PRS_3_Q20)
#1      2      3 
#40561 121679  40561 

#人年
fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male[male$PRO_PRS_3_Q20==1,],scale = 1)
fit2$pyears#
fit2$event#
#706/435555

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male[male$PRO_PRS_3_Q20==2,],scale = 1)
fit2$pyears#
fit2$event# 
#5108/1293931

fit2 <- pyears(Surv(PRC_difftime_new,PRC_total.y)~1,data=male[male$PRO_PRS_3_Q20==3,],scale = 1)
fit2$pyears#
fit2$event#
#3771/421144


#HR
model11 <- coxph(Surv(PRC_difftime_new,PRC_total.y)~as.factor(PRO_PRS_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
res=summary(model11)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                          coef exp(coef) se(coef)      z Pr(>|z|)
#as.factor(PRO_PRS_3_Q20)2            0.90430      2.47  0.04017 22.513 3.1e-112
#as.factor(PRO_PRS_3_Q20)3            1.75931      5.81  0.04112 42.786  0.0e+00

#                                       exp(coef) exp(-coef) lower .95 upper .95
#as.factor(PRO_PRS_3_Q20)2                2.47       0.40      2.28      2.67
#as.factor(PRO_PRS_3_Q20)3                5.81       0.17      5.36      6.30


#p for trend
model12 <- coxph(Surv(PRC_difftime_new,PRC_total.y)~as.numeric(PRO_PRS_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=male)
res=summary(model12)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                        coef exp(coef) se(coef)      z Pr(>|z|)
#as.numeric(PRO_PRS_3_Q20)            0.86828      2.38  0.01690 51.377  0.0e+00

#                                       exp(coef) exp(-coef) lower .95 upper .95
#as.numeric(PRO_PRS_3_Q20)                2.38       0.42      2.31      2.46



#乳腺癌  BC_PRS#
#female#
Q20=quantile(female$BC_PRS,seq(0.05,1,0.05))[[4]]
Q80=quantile(female$BC_PRS,seq(0.05,1,0.05))[[16]]

####对原始数据的PRS三分位，<=Q20，1；>=Q80，3
female$BC_PRS_3_Q20=0
female[female$BC_PRS<=Q20,]$BC_PRS_3_Q20=1
female[female$BC_PRS>Q20&female$BC_PRS<Q80,]$BC_PRS_3_Q20=2
female[female$BC_PRS>=Q80,]$BC_PRS_3_Q20=3
table(female$BC_PRS_3_Q20)
#1      2      3 
#47920 143758  47920 

#人年
fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female[female$BC_PRS_3_Q20==1,],scale = 1)
fit2$pyears#
fit2$event#
#1008/521612

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female[female$BC_PRS_3_Q20==2,],scale = 1)
fit2$pyears#
fit2$event#
#4758/1553572

fit2 <- pyears(Surv(BC_difftime_new,BC_total)~1,data=female[female$BC_PRS_3_Q20==3,],scale = 1)
fit2$pyears#
fit2$event#
#2430/511706


#HR
model11 <- coxph(Surv(BC_difftime_new,BC_total)~as.factor(BC_PRS_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
summary(model11)
res=summary(model11)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                       coef exp(coef) se(coef)     z Pr(>|z|)
#as.factor(BC_PRS_3_Q20)2             0.47273      1.60  0.03469 13.63  2.8e-42
#as.factor(BC_PRS_3_Q20)3             0.92901      2.53  0.03757 24.73 5.7e-135

#                                      exp(coef) exp(-coef) lower .95 upper .95
#as.factor(BC_PRS_3_Q20)2                 1.60       0.62      1.50       1.72
#as.factor(BC_PRS_3_Q20)3                 2.53       0.39      2.35       2.73

#p for trend
model12 <- coxph(Surv(BC_difftime_new,BC_total)~as.numeric(BC_PRS_3_Q20)+Age+Height_imp+Townsend_deprivation_index_imp_cat+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=female)
res=summary(model12,digit=10)
print(res$coefficients,digit=2)
print(res$conf.int,digit=2)
#                                       coef exp(coef) se(coef)     z Pr(>|z|)
#as.numeric(BC_PRS_3_Q20)             0.46239      1.59  0.01773 26.08 5.7e-150

#                                       exp(coef) exp(-coef) lower .95 upper .95
#as.numeric(BC_PRS_3_Q20)                 1.59       0.63      1.53       1.6

