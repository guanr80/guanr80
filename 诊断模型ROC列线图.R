#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af

#引用包
library(dplyr)
library(pROC)
library(ggplot2)
library(survival)
#BiocManager::install("regplot")
library(regplot)
#BiocManager::install("rms")
library(rms)
library(ggsci)
#BiocManager::install("survminer")
library(survminer)
#BiocManager::install("timeROC")
library(timeROC)
#BiocManager::install("ggDCA")
library(ggDCA)
library(limma)

inputFile="GSE2871.txt"       #表达矩阵
hub="hub.txt"        #核心基因

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
#rownames(rt)=rt[,1]
#exp=rt[,2:ncol(rt)]
exp=rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(rt)
sample=read.table("sample.txt",sep="\t",header=F,check.names=F)
colnames(sample)=c("ID","Type")
data=data[sample$ID,]
aSAH1=data[,read.table(hub, header=F, sep="\t", check.names=F)[,1]]
aSAH=cbind(sample,aSAH1)
#简单看一下ROC曲线AUC的情况
aflist=roc(Type~LOC103691092+Npw, data = aSAH)
g3 <- ggroc(aflist, size = 1.2,alpha=.6,)
g5=g3+ggsci::scale_color_lancet()
print(g5)
#诺曼图，高度比：8比10
############逻辑回归模型
dd <- datadist(aSAH)
options(datadist="dd")
fit <- lrm(formula = Type~LOC103691092+Npw, data =aSAH)
print(fit)
coef=as.data.frame(fit$coefficients)[-1,,drop=F]
coefout=cbind(ID=rownames(coef),coef)
write.table(coefout,file="coefficients.txt",sep="\t",quote=F,row.names = F)
#绘图
pdf(file="nomogram.pdf", width=9, height=7.5)
plot(nomogram(fit,fun.at = seq(0.05,0.95,0.05)),funlabel = "nomogram model")
dev.off()

plot(regplot(fit,plots=c("density","boxes"), observation=T, title="Prediction Nomogram", clickable=F, points=TRUE,droplines=TRUE))

nomoscore=predict(fit, data=t(aSAH))
aSAH$nomoscore=nomoscore
write.table(aSAH,file="nomoscore.txt",sep="\t",quote=F,row.names = F)

#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af
