#不同CPRS等级中有无多肿瘤家族史的比例（新构建CPRS）

#加载数据
load("imputed_final_20240120_UKB_database_for_analysis_202801_without_prevalent_cancer_male.rdata")
load("imputed_final_20240120_UKB_database_for_analysis_239598_without_prevalent_cancer_female.rdata")



#加载R包
library("ggplot2")
library("openxlsx")
library("data.table")
library("tidyverse")
library("writexl")
library("scales")
library("eoffice") 



###male###

#CPRS分组
Q5=quantile(male$CPRS,seq(0.05,1,0.05))[[1]]
Q10=quantile(male$CPRS,seq(0.05,1,0.05))[[2]]
Q20=quantile(male$CPRS,seq(0.05,1,0.05))[[4]]
Q25=quantile(male$CPRS,seq(0.05,1,0.05))[[5]]
Q75=quantile(male$CPRS,seq(0.05,1,0.05))[[15]]
Q80=quantile(male$CPRS,seq(0.05,1,0.05))[[16]]
Q90=quantile(male$CPRS,seq(0.05,1,0.05))[[18]]
Q95=quantile(male$CPRS,seq(0.05,1,0.05))[[19]]
Q30=quantile(male$CPRS,seq(0.05,1,0.05))[[6]]
Q40=quantile(male$CPRS,seq(0.05,1,0.05))[[8]]
Q50=quantile(male$CPRS,seq(0.05,1,0.05))[[10]]
Q60=quantile(male$CPRS,seq(0.05,1,0.05))[[12]]
Q70=quantile(male$CPRS,seq(0.05,1,0.05))[[14]]

####对CPRS五分位
male$CPRS_5=0
male[male$CPRS<=Q20,]$CPRS_5=1
male[male$CPRS>Q20&male$CPRS<=Q40,]$CPRS_5=2
male[male$CPRS>Q40&male$CPRS<=Q60,]$CPRS_5=3
male[male$CPRS>Q60&male$CPRS<=Q80,]$CPRS_5=4
male[male$CPRS>Q80,]$CPRS_5=5

table(male$CPRS_5)
#1     2     3     4     5 
#40561 40560 40560 40560 40560


#CPRS
table(male[male$CPRS_5==1,]$FH_total_2)   #30414 10147  0.2501664
table(male[male$CPRS_5==2,]$FH_total_2)   #29722 10838  0.2672091
table(male[male$CPRS_5==3,]$FH_total_2)   #29366 11194  0.2759862
table(male[male$CPRS_5==4,]$FH_total_2)   #28901 11659  0.2874507
table(male[male$CPRS_5==5,]$FH_total_2)   #28121 12439  0.3066815

dat=read.xlsx("D://2022.11.10   PRS & 家族史/20240120最新PRS/分析过程与结果/Sfigure 3/中间文件/CPRS_FH_male.xlsx")

p1 <- ggplot(dat,aes(x = CPRS,y = proportion))+
  geom_bar(stat = 'identity',aes(fill = FH),width = 0.9,position = position_stack(reverse = TRUE)) +
  theme(text=element_text(size=12,face = "bold"), #设置文字的字体字号（设置字体是为了确保汉字可以显示，字号和加粗请随意）
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=12))+
  scale_fill_manual(values=c("#00468B99","#AD002A99") ) +# 设置X轴文字大小
  scale_y_continuous(labels=percent)+# 使纵坐标呈现百分比
  theme(panel.grid = element_blank())
  
p1

ggsave("Proportion_of_multiple_FH_in_different_class_of_CPRS_of_male.pdf",print(p1),width = 8, height = 6, dpi = 600)
#551×380
topptx(p1,"Proportion_of_multiple_FH_in_different_class_of_CPRS_of_male.pptx",width = 5,height = 4)
#sgv width1000 height800



###female###

#CPRS分组
Q5=quantile(female$CPRS,seq(0.05,1,0.05))[[1]]
Q10=quantile(female$CPRS,seq(0.05,1,0.05))[[2]]
Q20=quantile(female$CPRS,seq(0.05,1,0.05))[[4]]
Q25=quantile(female$CPRS,seq(0.05,1,0.05))[[5]]
Q75=quantile(female$CPRS,seq(0.05,1,0.05))[[15]]
Q80=quantile(female$CPRS,seq(0.05,1,0.05))[[16]]
Q90=quantile(female$CPRS,seq(0.05,1,0.05))[[18]]
Q95=quantile(female$CPRS,seq(0.05,1,0.05))[[19]]
Q30=quantile(female$CPRS,seq(0.05,1,0.05))[[6]]
Q40=quantile(female$CPRS,seq(0.05,1,0.05))[[8]]
Q50=quantile(female$CPRS,seq(0.05,1,0.05))[[10]]
Q60=quantile(female$CPRS,seq(0.05,1,0.05))[[12]]
Q70=quantile(female$CPRS,seq(0.05,1,0.05))[[14]]

####对原始数据的CPRS五分位
female$CPRS_5=0
female[female$CPRS<=Q20,]$CPRS_5=1
female[female$CPRS>Q20&female$CPRS<=Q40,]$CPRS_5=2
female[female$CPRS>Q40&female$CPRS<=Q60,]$CPRS_5=3
female[female$CPRS>Q60&female$CPRS<=Q80,]$CPRS_5=4
female[female$CPRS>Q80,]$CPRS_5=5

table(female$CPRS_5)
#1     2     3     4     5 
#47920 47919 47920 47919 47920


#CPRS
table(female[female$CPRS_5==1,]$FH_total_2)   #35545 12375  0.2582429
table(female[female$CPRS_5==2,]$FH_total_2)   #33853 14066  0.293537
table(female[female$CPRS_5==3,]$FH_total_2)   #33305 14615  0.3049875
table(female[female$CPRS_5==4,]$FH_total_2)   #32778 15141  0.3159707
table(female[female$CPRS_5==5,]$FH_total_2)   #32053 15867  0.3311144

dat=read.xlsx("D://2022.11.10   PRS & 家族史/20240120最新PRS/分析过程与结果/Sfigure 3/中间文件/CPRS_FH_female.xlsx")

p2 <- ggplot(dat,aes(x = CPRS,y = proportion))+
  geom_bar(stat = 'identity',aes(fill = FH),width = 0.9,position = position_stack(reverse = TRUE)) +
  theme(text=element_text(size=12,face = "bold"), #设置文字的字体字号（设置字体是为了确保汉字可以显示，字号和加粗请随意）
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=12))+
  scale_fill_manual(values=c("#00468B99","#AD002A99") ) +# 设置X轴文字大小
  scale_y_continuous(labels=percent)+# 使纵坐标呈现百分比
  theme(panel.grid = element_blank())
p2

ggsave("Proportion_of_multiple_FH_in_different_class_of_CPRS_of_female.pdf",print(p2),width = 8, height = 6, dpi = 600)
topptx(p2,"Proportion_of_multiple_FH_in_different_class_of_CPRS_of_female.pptx",width = 5,height = 4)
