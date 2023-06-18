#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af


#bioclite()换个镜像
#BiocManager::install("randomcoloR")
library(randomcoloR)
library(venn) 

#文件名
A="WGCNA"
B="DIFF"

#构建一个列表
geneList=list()
rt=read.table(paste0(A,".txt"),header=F,sep="\t",check.names=F)      
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)      
uniqGene=unique(geneNames)                 
geneList[[A]]=uniqGene                   
uniqLength=length(uniqGene)
print(paste("1",uniqLength,sep=" "))
rt=read.table(paste0(B,".txt"),header=F,sep="\t",check.names=F)    
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)       
uniqGene=unique(geneNames)                
geneList[[B]]=uniqGene
uniqLength=length(uniqGene)
print(paste("3",uniqLength,sep=" "))

mycol <- distinctColorPalette(100)

pdf(file="disease.pdf",width=10,height=10)  

venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F)
attach(mtcars)

plot(wt, mpg)
dev.off()
dev.off()

intersectGenes=Reduce(intersect,geneList)
intersectGenes
write.table(file="disease.txt",intersectGenes,sep="\t",quote=F,col.names=F,row.names=F) 

#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af
