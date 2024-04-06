#有无家族史人群PRS分布

#######画PRS分布图  选择用平滑曲线

#加载数据
load("imputed_final_20240120_UKB_database_for_analysis_202801_without_prevalent_cancer_male.rdata")
load("imputed_final_20240120_UKB_database_for_analysis_239598_without_prevalent_cancer_female.rdata")


#加载R包
library("data.table")
library("ggsci")
library("ggplot2")
library("eoffice")
library("Cairo")


####肠癌 COL_PRS
#男性
#肠癌家族史二分类0，1，   FH_Bowel_2

table(male$FH_Bowel_2)
#0      1 
#180048  22753 

dat=male

dat$colour=0
dat$colour[dat$FH_Bowel_2==0]="#00468B99"
dat$colour[dat$FH_Bowel_2==1]="#AD002A99"


p2 <- ggplot(dat)+geom_density(aes(x = COL_PRS, fill = colour),alpha=0.35)+
  xlab("Polygenic Risk Score")+ylab("Proportion")+
  scale_fill_lancet(labels = c("No family cancer","With family cancer"))+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),legend.position=c(0.2,0.85),legend.title=element_blank(),
        legend.text = element_text(size = 16),
        legend.background = element_blank())
p2


ggsave("Distribution_of_male_colorectal_prs_with_or_without_family history.pdf", plot = p2, width = 8, height = 6, dpi = 600)


#女性#
dat=female

dat$colour=0
dat$colour[dat$FH_Bowel_2==0]="#00468B99"
dat$colour[dat$FH_Bowel_2==1]="#AD002A99"


###平滑曲线
p2 <- ggplot(dat)+geom_density(aes(x = COL_PRS, fill = colour),alpha=0.35)+
  xlab("Polygenic Risk Score")+ylab("Proportion")+
  scale_fill_lancet(labels = c("No family cancer","With family cancer"))+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),legend.position=c(0.2,0.85),legend.title=element_blank(),
        legend.text = element_text(size = 16),
        legend.background = element_blank())

p2
ggsave("Distribution_of_female_colorectal_prs_with_or_without_family history.pdf", plot = p2, width = 8, height = 6, dpi = 600)






#肺癌
#male#
dat=male

dat$colour=0
dat$colour[dat$FH_Lung_2==0]="#00468B99"
dat$colour[dat$FH_Lung_2==1]="#AD002A99"


###平滑曲线
p2 <- ggplot(dat)+geom_density(aes(x = PRS_LC, fill = colour),alpha=0.35)+
  xlab("Polygenic Risk Score")+ylab("Proportion")+
  scale_fill_lancet(labels = c("No family cancer","With family cancer"))+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),legend.position=c(0.18,0.85),legend.title=element_blank(),
        legend.text = element_text(size = 16),
        legend.background = element_blank())
p2


ggsave("Distribution_of_male_lung_prs_with_or_without_family history.pdf", plot = p2, width = 8, height = 6, dpi = 600)



#female#
dat=female

dat$colour=0
dat$colour[dat$FH_Lung_2==0]="#00468B99"
dat$colour[dat$FH_Lung_2==1]="#AD002A99"



###平滑曲线
p2 <- ggplot(dat)+geom_density(aes(x = PRS_LC, fill = colour),alpha=0.35)+
  xlab("Polygenic Risk Score")+ylab("Proportion")+
  scale_fill_lancet(labels = c("No family cancer","With family cancer"))+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),legend.position=c(0.17,0.85),legend.title=element_blank(),
        legend.text = element_text(size =16),
        legend.background = element_blank())  # 移除图例的背景
p2



ggsave("Distribution_of_female_lung_prs_with_or_without_family history.pdf", plot = p2, width = 8, height = 6, dpi = 600)



#前列腺癌#
dat=male

dat$colour=0
dat$colour[dat$FH_Prostate_2==0]="#00468B99"
dat$colour[dat$FH_Prostate_2==1]="#AD002A99"



###平滑曲线
p2 <- ggplot(dat)+geom_density(aes(x = PRS_PRC, fill = colour),alpha=0.35)+
  xlab("Polygenic Risk Score")+ylab("Proportion")+
  scale_fill_lancet(labels = c("No family cancer","With family cancer"))+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),legend.position=c(0.2,0.85),legend.title=element_blank(),
        legend.text = element_text(size =16),
        legend.background = element_blank())  # 移除图例的背景
p2


ggsave("Distribution_of_male_prostate_prs_with_or_without_family history.pdf", plot = p2, width = 8, height = 6, dpi = 600)



#乳腺癌#
dat=female

dat$colour=0
dat$colour[dat$FH_Breast_2==0]="#00468B99"
dat$colour[dat$FH_Breast_2==1]="#AD002A99"



###平滑曲线
p2 <- ggplot(dat)+geom_density(aes(x = PRS_BRC, fill = colour),alpha=0.35)+
  xlab("Polygenic Risk Score")+ylab("Proportion")+
  scale_fill_lancet(labels = c("No family cancer","With family cancer"))+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),legend.position=c(0.2,0.85),legend.title=element_blank(),
        legend.text = element_text(size =16),
        legend.background = element_blank())  # 移除图例的背景

p2


ggsave("Distribution_of_female_breast_prs_with_or_without_family history.pdf", plot = p2, width = 8, height = 6, dpi = 600)


