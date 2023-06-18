#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af

##############################################################5.ppi
library(ggsci) 
library(tidyverse) 
library(clusterProfiler) 
library(org.Hs.eg.db)
library(org.Mm.eg.db)
#BiocManager::install("STRINGdb")
library(STRINGdb) 

library(igraph) 
library(ggraph) 
library(tidygraph) 
library(cluster) 
library(factoextra) 
mycol <- pal_aaas("default", alpha = 0.8)(10) #ggsci颜色选取

#读入
gene=read.table("disease.txt",header=F,sep="\t",check.names=F)[,1]

#### 1.STRING PPI ####
#string_db <- STRINGdb$new(version="11", species=9606, #9606为智人，10090为小鼠
                          #input_directory="",score_threshold=400) #400为置信度
string_db <- STRINGdb$new(version="11", species=10090, #9606为智人，10090为小鼠
                          input_directory="",score_threshold=400) #400为置信度
#gene <- gene %>% bitr(fromType = "SYMBOL", 
                      #toType = "ENTREZID", 
                      #OrgDb = "org.Hs.eg.db", #小鼠将Hs改为Mm
                      #drop = T)
gene <- gene %>% bitr(fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = "org.Mm.eg.db", #小鼠将Hs改为Mm
                      drop = T)
data_mapped <- gene %>% string_db$map(my_data_frame_id_col_names = "ENTREZID", 
                                      removeUnmappedRows=T,quiet=T) #较耗时

pdf(file=paste0("6.","_STRING.pdf"),height=10,width=10) #保存图片
string_db$plot_network(data_mapped$STRING_id) 
dev.off()
hit<-data_mapped$STRING_id
info<-string_db$get_interactions(hit)
links <- info %>% 
  mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
  mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
  dplyr::select(from, to , last_col()) %>% 
  dplyr::rename(weight = combined_score) %>% 
  distinct(from,to,weight)
write.table(links,paste0("6.","_links.txt"),sep="\t",quote=F,row.names=F) #保存结果


#### 2.GGRAPH PPI ####
links=read.table(paste0("6.","_links.txt"),sep="\t",header=T) #输入结果
# st=dir()[grepl("string_interactions_short.tsv",dir())]
# links=read.table(st)[,c(1,2,13)]
# links=links %>% rename(from=V1,to=V2,weight=V13)
# links$weight=links$weight*1000
#从STRING网站下载的表格放在同文件夹下
links_2 <- links %>% mutate(from_c = count(., from)$n[match(from, count(., from)$from)]) %>%
  mutate(to_c = count(., to)$n[match(to, count(., to)$to)]) %>%
  filter(!(from_c <= 2 & to_c <= 2)) %>% #去除游离节点
  dplyr::select(1,2,3) 
nodes_2 <- links_2 %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
net_2 <- igraph::graph_from_data_frame(d=links_2,vertices=nodes_2,directed = F)
igraph::V(net_2)$deg <- igraph::betweenness(net_2) #绘制基于betweenness值筛选hub基因PPI
igraph::V(net_2)$size <- igraph::betweenness(net_2)
igraph::V(net_2)$hub <- ifelse(betweenness(net_2)>mean(betweenness(net_2)),
                               pal_locuszoom("default",alpha=0.8)(7)[1],
                               pal_locuszoom("default",alpha=0.8)(7)[3])
igraph::E(net_2)$width <- igraph::E(net_2)$weight
pdf(paste0("6.","_PPI_hub.pdf"),height=6,width=8)
ggraph(net_2,layout = 'stress')+
  geom_edge_link(color="grey",aes(edge_width=width),alpha=0.4,show.legend=F)+
  scale_edge_width(guide = "none",range=c(0.6,1.2))+ #调节基于置信度的线条粗细的最大最小值
  geom_node_point(aes(size=2*log((V(net_2)$size+mean(V(net_2)$size))/mean(V(net_2)$size)+1)), #看到2*了么，可修改节点大小
                  show.legend=F,color=V(net_2)$hub, alpha=0.8)+
  geom_node_text(aes(label=name,size=sqrt(log((V(net_2)$size+mean(V(net_2)$size))/mean(V(net_2)$size)+1))),
                 show.legend=F,repel=F)+
  theme_graph()
dev.off()

#保存蛋白互作网络的核心基因(选)
aa=as.data.frame(betweenness(net_2))
aa=aa[aa$`betweenness(net_2)`>mean(aa$`betweenness(net_2)`),,drop=F]
aa=rownames(aa)
write.table(aa,file="PPI.hub.txt",sep="\t",quote=F,col.names = F,row.names = F)
