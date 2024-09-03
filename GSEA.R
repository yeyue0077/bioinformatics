library(ggplot2)
#BiocManager::install("limma")
library(limma)
library(pheatmap)
library(ggsci)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(readxl)

#输入文件
data_all<-read.csv("合并数据集表达谱-2.csv",row.names = 1)
colnames(data_all)<-gsub("\\.","-",colnames(data_all))
Group<- read.csv("合并数据集分组-2.csv")
colnames(Group)
isr_ids <- Group$ID[Group$group == "ISR"]  
data <-data_all[,isr_ids]


data <- as.matrix(data)
gene <-read_xlsx("gene.xlsx")
colnames(gene)[1]<-"GENE"


sgene= gene$GENE[2]      #这里输入进行单基因GSEA的基因名称
#下面我们对表达矩阵处理一下，先读入
sgene
data=avereps(data)


#按基因表达分组
group <- ifelse(data[c(sgene),]>median(data[c(sgene),]), "High", "Low")   
group <- factor(group,levels = c("High","Low"))
#我们这里是是没有做差异分析的矩阵，我们需要找到与我们想要基因相关的，并且我们的基因也需要有差异才可以真实性
#接下来：
#差异分析
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(High-Low,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
deg=topTable(fit2,adjust='fdr',number=nrow(data))
Diff=deg
#保存单基因分组的所有基因差异结果
DIFFOUT=rbind(id=colnames(Diff),Diff)
write.table(DIFFOUT,file=paste0("ANKRD13A.","DIFF_all.xls"),sep="\t",quote=F,col.names=F)


#展示差异最大的前30个基因
Diff=Diff[order(as.numeric(as.vector(Diff$logFC))),]
diffGene=as.vector(rownames(Diff))
diffLength=length(diffGene)
afGene=c()
if(diffLength>(60)){
  afGene=diffGene[c(1:30,(diffLength-30+1):diffLength)]
}else{
  afGene=diffGene
}
afExp=data[afGene,]
#我们对分组进行标记，
#分组标签
Type1=as.data.frame(group)
Type1=Type1[order(Type1$group,decreasing = T),,drop=F]
Type=Type1[,1]
names(Type)=rownames(Type1)
Type=as.data.frame(Type)
#分组标签的注释颜色
anncolor=list(Type=c(High="#D77071",Low="#6888F5"  ))

pdf(file=paste0("LRRK2.", "DIFF_heatmap.pdf"),height=7,width=6)
pheatmap(afExp[,rownames(Type1)],                                                                      #热图数据
         annotation=Type,                                                            #分组
         color = colorRampPalette(c("#6888F5","white","#D77071"))(50),     #热图颜色
         cluster_cols =F,                                                           #不添加列聚类树
         show_colnames = F,                                                         #展示列名
         scale="row", 
         fontsize = 10,
         fontsize_row=6,
         fontsize_col=8,
         annotation_colors=anncolor
)
dev.off()


#火山图差异标准设置，这个小伙伴自行设置，并不是一定要绘制
# 设置阈值  
adjP <- 0.05  
aflogFC <- 0.5  

# 创建一个新的变量来指定颜色  
Diff$Significant <- ifelse(Diff$P.Value < adjP & Diff$logFC > aflogFC, "Up",  
                           ifelse(Diff$P.Value < adjP & Diff$logFC < -aflogFC, "Down", "NotSig"))  

# 开始绘制火山图  
p <- ggplot(Diff, aes(logFC, -log10(P.Value))) +  
  geom_point(aes(col = Significant), size = 4) +               # 根据新的color变量设置点的颜色  
  scale_color_manual(values = c("Up" = "#D77071", "Down" = "#6888F5", "NotSig" = "grey")) + # 设置颜色  
  labs(title = "Volcano Plot") +  # 添加标题  
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) +  
  geom_hline(yintercept = -log10(adjP), colour = "gray", linetype = "twodash", size = 1) +  
  geom_vline(xintercept = aflogFC, colour = "gray", linetype = "twodash", size = 1) +  
  geom_vline(xintercept = -aflogFC, colour = "gray", linetype = "twodash", size = 1)  

# 查看火山图  
p

#添加基因点标记，按照,可自行根据差异分析的结果进行标记
point.Pvalue=0.01
point.logFc=1.5
deg$symbol=rownames(deg)
highlighted_genes <- deg[deg$P.Value < point.Pvalue & abs(deg$logFC) > point.logFc, ]  
# 开始绘制火山图  
p <- ggplot(Diff, aes(logFC, -log10(P.Value))) +  
  geom_point(aes(col = Significant), size = 4) +               # 根据新的color变量设置点的颜色  
  scale_color_manual(values = c("Up" = "#D77071", "Down" = "#6888F5", "NotSig" = "grey")) + # 设置颜色  
  labs(title = "Volcano Plot") +  # 添加标题  
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) +  
  geom_hline(yintercept = -log10(adjP), colour = "gray", linetype = "twodash", size = 1) +  
  geom_vline(xintercept = aflogFC, colour = "gray", linetype = "twodash", size = 1) +  
  geom_vline(xintercept = -aflogFC, colour = "gray", linetype = "twodash", size = 1) +
  
  # 添加满足条件的基因的标签  
  geom_text(data = highlighted_genes, aes(label = symbol), # 使用满足条件的基因数据框  
            vjust = -0.5, # 垂直调整文本位置  
            color = "black") # 设置文本颜色，可以根据需要调整  

# 查看火山图  
p



########转换ID#########
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort <- DEG %>% 
  arrange(desc(logFC))

geneList = data_all_sort$logFC #把foldchange按照从大到小提取出来
names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(geneList)

#开始富集分析
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 10,
               maxGSSize    = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none" )
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
# 获取result中按enrichmentScore排序后的前几个基因集的行名（通常是通路的名称）  
# 注意：这里只是获取了行名，并没有实际地修改或选择数据
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af=as.data.frame(kk2@result)
write.table(af,file=paste0("LRRK2.","all_GSEA.xls"),sep="\t",quote=F,col.names=T)#这里我们把结果保存一下


#排序后分别取GSEA结果的前5个和后5个
num=5
pdf(paste0("LRRK2.","down_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)])
dev.off()
pdf(paste0("LRRK2.","up_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore),num)])
dev.off()

#排序后取前5个和后5个一起展示
num=5
pdf(paste0("LRRK2.","all_GSEA.pdf"),width = 10,height = 10)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))])
dev.off()
#NES绝对值最大的5个展示#
num=5
pdf(paste0("LRRK2.","first_abs_NES_5_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[order(abs(kk2@result$NES), decreasing = TRUE)][1:5])
dev.off()

##某一通路单独###
gseaplot2(kk2,
          title = "hsa05134",  #设置标题
          "hsa05134", #绘制hsa05134通路的结果，通路名称与编号对应
          color="red", #线条颜色
          base_size = 20, #基础字体的大小
          subplots = 1:3, 
          pvalue_table = T) # 显示p值
