######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")


#???ð?
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
expFile="ssgseaOut.txt"      #?????????ļ?
riskFile="sample.txt"       #?????????ļ?
geneFile="immuneset.txt"       #???߼??????Ļ????ļ?
setwd("D:\\biowolf\\FerrLnc\\29.checkpoint")     #???ù???Ŀ¼

#??ȡ?????????ļ?,???????ݽ??д???
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
	
#??ȡ?????ļ?
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data),as.vector(gene[,1]))
data=t(data[sameGene,])

#?ϲ?????
risk=read.table(riskFile, sep="\t", header=F, check.names=F, row.names=1)
colnames(risk)="data"
sameSample=intersect(row.names(data),row.names(risk))
rt1=cbind(data[sameSample,],risk[sameSample,])
colnames(rt1)[ncol(rt1)]="data"
rt1=rt1[,c(sameGene,"data")]
rt1=as.data.frame(rt1)
rt1[,1:(ncol(rt1)-1)]=lapply(rt1[,1:(ncol(rt1)-1)],as.numeric)
#??ȡ?????????Ļ???
sigGene=c()
for(i in colnames(rt1)[1:(ncol(rt1)-1)]){
	if(sd(rt1[,i])<0.001){next}
	wilcoxTest=wilcox.test(rt1[,i] ~ rt1[,"data"])
	pvalue=wilcoxTest$p.value
	if(wilcoxTest$p.value<0.05){
		sigGene=c(sigGene, i)
	}
}
sigGene=c(sigGene, "data")
rt1=rt1[,sigGene]

#??????ת????ggplot2?????ļ?
rt1=melt(rt1,id.vars=c("data"))
colnames(rt1)=c("data","Gene","Expression")
	
#???ñȽ???
group=levels(factor(rt1$data))
rt1$data=factor(rt1$data, levels=c("MI","Healthy"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}
	
#????????ͼ
boxplot=ggboxplot(rt1, x="Gene", y="Expression", color ="data",
				  xlab="",
				  ylab="Score",
				  legend.title="Cluster",
				  add = "jitter",
				  width=0.8,
				  palette = ggsci::pal_npg()(2) )+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=data),
	method="wilcox.test",
	symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 0.2,1), symbols=c("***", "**", "*", "#","ns")), label="p.signif")
	
#????ͼƬ
pdf(file="MI.diff.pdf",  width=10, height=8)
print(boxplot)
dev.off()


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056
