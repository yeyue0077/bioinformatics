#8.ssGSEA##################################################################################################
###a.读取data_all和Group##############################################################
library(data.table)
library(dplyr)
library(GSVA)
#BiocManager::install("GSVA")
library(ComplexHeatmap)
library(ggplot2)
library(corrgram)
library(ggthemes)
library(ggpubr)
library(readxl)
#读取
data_all<-read.csv("合并数据集表达谱-2.csv",row.names = 1)
colnames(data_all)<-gsub("\\.","-",colnames(data_all))
# row.names(data_all) <- data_all[,1]
# data_all <- data_all[,-1]

Group<- read.csv("合并数据集分组-2.csv")
colnames(Group)
colnames(Group)<-c("ID","group")
# colnames(Group)<-c("ID","group")
# colnames(tcga)<-gsub("\\.","-",colnames(tcga))#把列名中的点换成"-"
data_all<-data_all[,Group$ID]# 按照分组信息筛选表达数据

dir.create("output")
#####end#################################################################################

###b.获取细胞浸润丰度####################################################################
data<-data_all# 将表达数据赋值给变量data
group<-Group# 将分组信息赋值给变量group
cell_gene<-read.table("Immune marker genes.txt",sep="\t")# 读取免疫细胞标记基因
cell<-unique(as.character(cell_gene[,2]))# 获取所有免疫细胞类型
cell_gene_type=list()# 创建一个空列表用于存储每种细胞类型的标记基因
# 将标记基因按照细胞类型分组
for (i in 1:length(cell)){
  b=list(as.character(cell_gene[(as.character(cell_gene$V2))==cell[i],1]))
  cell_gene_type[i]=(b)
  names(cell_gene_type)[i]=cell[i]
}
# 使用ssGSEA算法计算免疫细胞浸润分数
gsva_es1<-gsva(as.matrix(data),cell_gene_type,method="ssgsea",abs.ranking=F)
#############处理警告信息，新的调用方式####
params <- ssgseaParam(as.matrix(data),cell_gene_type)
gsva_es1 <- gsva(params)
############################################
write.csv(gsva_es1, "output/1-immu_ssgsea.csv")
#####end#################################################################################

###c.绘制免疫细胞在不同分组中的分组比较图#########################################################
cb2<-as.data.frame(t(gsva_es1)) # 转置浸润分数矩阵
# cb2<-cb2[,which(colnames(cb2)%in% rownames(output))]
immu<- rep(colnames(cb2),each=nrow(cb2))  # 创建免疫细胞变量
immu <- factor(immu) #组别因子化
a<-c(group$group)# 获取分组信息
group <- rep(a,ncol(cb2))  # 重复分组信息以匹配免疫细胞数量
group <- factor(group,levels = c('Control','ISR')) # 将分组变量转换为因子并指定顺序
value <- c() # 创建空向量用于存储浸润分数
# 将浸润分数转换为长格式
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value) #生成数据框

#免疫细胞的分组比较图
immu_pvalue <- compare_means(value~group,data = Data,group.by = "immu_cell")
immu_pvalue <- immu_pvalue %>% arrange(immu_cell)
write.csv(immu_pvalue, "output/boxplot-immu_pvalue.csv")

p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA,linewidth=0.1,fatten=1,position=position_dodge(0.85))+
  theme_bw()+scale_fill_manual(values = c("#6888F5","#D77071"))+#颜色顺序同前面group因子化顺序一致。
  theme(panel.grid = element_blank())+
  stat_compare_means(symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")),label = "p.signif",size=1.5,hjust = 0.5,vjust = 0)+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('Infiltration Abundance')+labs(fill='Group')+
  theme(axis.text = element_text(size=6,colour = "black"),
        legend.text = element_text(size=6,colour = "black"),
        legend.title = element_text(size=6,colour = "black"),
        axis.title = element_text(size=6,colour = "black"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),legend.box.spacing = unit(1,"pt"),
        legend.margin = ggplot2::margin(1,0,0,0),
        plot.margin = ggplot2::margin(0,2,0,1)
  )
p
ggsave("output/3-immu_barplot2.pdf",plot=p,width = 16,height = 7,units = "cm")
dev.off()
#####end#################################################################################

########全样本，所有免疫细胞间相关性#######
cells <- immu_pvalue$immu_cell %>% as.character()# 获取所有免疫细胞名称
col1 <- colorRampPalette(colors = c("#6888F5","white","#D77071"),space="Lab")#定义颜色,相关性排序：低，中，高
immu<-t(gsva_es1)
write.csv(immu, "output/immu-pvalue.csv")
library(corrplot)
#method参数圆"circle", 方块"square", 椭圆"ellipse", 数字"number", 阴影"shade", 颜色"color", 饼形图"pie"
pdf("output/4-immu_corrplot.pdf",width = 10/2.54,height = 8/2.54,pointsize = 6)

corrplot(corr =cor(immu,method="spearman"),order="alphabet",type="upper",tl.pos="tp",method="color",#method参数换成上面几个中的一个就可以
         tl.col = "black",col=col1(50),cl.cex = 1,tl.cex = 1,cl.ratio = 0.2,mar = c(1,2,1,0))
corrplot(corr = cor(immu,method="spearman"),add=TRUE, type="lower", method="pie",order="alphabet",
         diag=FALSE,tl.pos="n", cl.pos="n",col=col1(50),mar = c(1,2,1,0))#"black"
dev.off()
corr =cor(immu,method="spearman")
write.csv(corr,"output/4_cor.csv")

###d.可视化ISR分组中的免疫细胞间（差异细胞）相关性#####################################################

immu_pvalue_sig <- filter(immu_pvalue,immu_pvalue$p.signif != "ns")
cells <- immu_pvalue_sig$immu_cell %>% as.character()
col1 <- colorRampPalette(colors = c("#6888F5","white","#D77071"),space="Lab")#定义颜色,相关性排序：低，中，高
immu<-t(gsva_es1[which((rownames(gsva_es1) %in% cells)),which(a == "ISR")])
write.csv(immu, "output/immu-pvalue0.05-ISR.csv")
library(corrplot)
#method参数圆"circle", 方块"square", 椭圆"ellipse", 数字"number", 阴影"shade", 颜色"color", 饼形图"pie"
pdf("output/4-immu_corrplot_ISR.pdf",width = 10/2.54,height = 8/2.54,pointsize = 6)

corrplot(corr =cor(immu,method="spearman"),order="alphabet",type="upper",tl.pos="tp",method="color",#method参数换成上面几个中的一个就可以
         tl.col = "black",col=col1(50),cl.cex = 1,tl.cex = 1,cl.ratio = 0.2,mar = c(1,2,1,0))
corrplot(corr = cor(immu,method="spearman"),add=TRUE, type="lower", method="pie",order="alphabet",
         diag=FALSE,tl.pos="n", cl.pos="n",col=col1(50),mar = c(1,2,1,0))#"black"
dev.off()
corr =cor(immu,method="spearman")
write.csv(corr,"output/4_cor_ISR.csv")
#####end#################################################################################

###e.可视化Control分组中的免疫细胞间相关性#####################################################
#取分组比较图中具有统计学差异的免疫细胞与High样本的相关性矩阵
immu<-t(gsva_es1[which((rownames(gsva_es1) %in% cells)),which(a == "Control")])
write.csv(immu, "output/immu-pvalue0.05-Control.csv")

#method参数圆"circle", 方块"square", 椭圆"ellipse", 数字"number", 阴影"shade", 颜色"color", 饼形图"pie"
pdf("output/5-immu_corrplot-Control.pdf",width = 10/2.54,height = 8/2.54,pointsize = 6)

corrplot(corr =cor(immu,method="spearman"),order="alphabet",type="upper",tl.pos="tp",method="color",#method参数换成上面几个中的一个就可以
         tl.col = "black",col=col1(50),cl.cex = 1,tl.cex = 1,cl.ratio = 0.2,mar = c(1,2,1,0))
corrplot(corr = cor(immu,method="spearman"),add=TRUE, type="lower", method="pie",order="alphabet",
         diag=FALSE,tl.pos="n", cl.pos="n",col=col1(50),mar = c(1,2,1,0))#"black"  
dev.off()
corr =cor(immu,method="spearman")
write.csv(corr,"output/5_cor_Control.csv")
#####end#################################################################################




group1<- Group
bb<-read_xlsx("gene.xlsx")
colnames(bb)[1]<-"GENE"
data1<-data[bb$GENE,]
library(ggplot2)
library(dplyr)
library(Hmisc)
library("Rmisc")
library("plyr")

###f.绘制 ISR组细胞（所有免疫细胞）与基因间相关性点图######################################################

# 获取所有免疫细胞名称
cells <- immu_pvalue$immu_cell %>% as.character()
# 构建数据矩阵：ISR组的免疫细胞浸润分数和基因表达数据
immu1<-cbind(t(gsva_es1[which((rownames(gsva_es1)%in% cells)),which(group1$group == "ISR")]),t(data1[,which(group1$group == "ISR")]))
# 初始化相关性结果数据框
cor<-data.frame()
# 循环计算相关性
for (mm in (nrow(immu_pvalue)+1):ncol(immu1)) {# 遍历每个基因（列）
  cor1<-data.frame(0,0,0,0) # 初始化临时数据框
  for (i in 1:nrow(immu_pvalue)) {# 遍历每个免疫细胞（行）
    c<-rcorr(immu1[,i],immu1[,mm],type ="spearman")# 计算Spearman相关性
    cor1[i,1]<-c$r[2]# 存储相关系数
    cor1[i,2]<-c$P[2] # 存储P值
    cor1[i,3]<-colnames(immu1)[i]# 存储免疫细胞名称
    cor1[i,4]<-colnames(immu1)[mm] # 存储基因名称
  }
  cor<-rbind(cor,cor1) # 将结果添加到cor数据框
}
colnames(cor)<-c("correlation","Pvalue","immu_cell","gene")

# 挑选出ISR样本中所有数据（注释中写的是Pvalue<0.05，但实际代码没有筛选）
cor4<-cor#[cor$Pvalue<0.05,]
cor4<-na.omit(cor4)# 删除缺失值
cor4<-cor4[order(cor4$correlation,decreasing = F),]# 按相关性系数升序排序
# 处理P值等于0的情况
CorPvalue <- cor4$Pvalue[which(cor4$Pvalue != 0)]# 获取非零P值
minCorPvalue <- min(CorPvalue)/2# 计算最小非零P值的一半
# 将P值等于0的替换为最小非零P值的一半
for (i in c(1:length(cor4$Pvalue))){
  if (cor4$Pvalue[i] == 0){
    cor4$Pvalue[i]=minCorPvalue
  }
}
write.csv(cor4, "output/6-immu_cor-cells-genes-ISR.csv")
y=factor(cor4$gene)
x=factor(cor4$immu_cell)
a<-(-log10(cor4$Pvalue))
p<-ggplot(cor4,aes(x,y))+ geom_point(aes(size=a,color=correlation),shape=20)+
  scale_colour_gradient(high="#D77071",low="#6888F5")+
  scale_size_continuous(range = c(0,3))+
  # theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Low")+
  theme_classic()+labs(size="-log10Pvalue",x="",y="")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(colour = "black",size=6),
        legend.text = element_text(colour = "black",size=6),
        legend.title = element_text(colour = "black",size=6),
        legend.key.size = unit(0.2,"cm"),legend.spacing = unit(0.1, "cm"),legend.margin = ggplot2::margin(6,2,0,0),
        legend.background = element_blank(),plot.background = element_blank(),
        plot.margin = ggplot2::margin(6,0,0,0))+
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 2))
p
ggsave(plot = p,file="output/6-immu-cell-ISR.pdf",width =8,height =6,units = "cm")
#####end#################################################################################


###f.绘制 ISR组细胞（差异细胞）与基因间相关性点图######################################################
cells <- immu_pvalue_sig$immu_cell %>% as.character()
immu1<-cbind(t(gsva_es1[which((rownames(gsva_es1)%in% cells)),which(group1$group == "ISR")]),t(data1[,which(group1$group == "ISR")]))
cor<-data.frame()
for (mm in (nrow(immu_pvalue_sig)+1):ncol(immu1)) {#目的基因所在列
  cor1<-data.frame(0,0,0,0)
  for (i in 1:nrow(immu_pvalue_sig)) {#免疫细胞所在列
    c<-rcorr(immu1[,i],immu1[,mm],type ="spearman")
    
    cor1[i,1]<-c$r[2]
    cor1[i,2]<-c$P[2]
    cor1[i,3]<-colnames(immu1)[i]
    cor1[i,4]<-colnames(immu1)[mm]
  }
  cor<-rbind(cor,cor1)
}
colnames(cor)<-c("correlation","Pvalue","immu_cell","gene")

#挑选出ISR样本中Pvalue<0.05的数据
cor4<-cor#[cor$Pvalue<0.05,]
cor4<-na.omit(cor4)
cor4<-cor4[order(cor4$correlation,decreasing = F),]

CorPvalue <- cor4$Pvalue[which(cor4$Pvalue != 0)]
minCorPvalue <- min(CorPvalue)/2
for (i in c(1:length(cor4$Pvalue))){
  if (cor4$Pvalue[i] == 0){
    cor4$Pvalue[i]=minCorPvalue
  }
}
write.csv(cor4, "output/6-immu_cor-cells-genes-ISR.csv")
y=factor(cor4$gene)
x=factor(cor4$immu_cell)
a<-(-log10(cor4$Pvalue))
p<-ggplot(cor4,aes(x,y))+ geom_point(aes(size=a,color=correlation),shape=20)+
  scale_colour_gradient(high="#FF6A6A",low="#DADDFC")+
  scale_size_continuous(range = c(0,3))+
  # theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Low")+
  theme_classic()+labs(size="-log10Pvalue",x="",y="")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(colour = "black",size=6),
        legend.text = element_text(colour = "black",size=6),
        legend.title = element_text(colour = "black",size=6),
        legend.key.size = unit(0.2,"cm"),legend.spacing = unit(0.1, "cm"),legend.margin = ggplot2::margin(6,2,0,0),
        legend.background = element_blank(),plot.background = element_blank(),
        plot.margin = ggplot2::margin(6,0,0,0))+
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 2))
p
ggsave(plot = p,file="output/6-immu-cell-ISR.pdf",width =8,height =6,units = "cm")
#####end#################################################################################


###g.绘制Control组细胞与基因间相关性点图######################################################
immu1<-cbind(t(gsva_es1[which((rownames(gsva_es1)%in% cells)),which(group1$group == "Control")]),t(data1[,which(group1$group == "Control")]))
cor<-data.frame()
for (mm in (nrow(immu_pvalue_sig)+1):ncol(immu1)) {#目的基因所在列
  cor1<-data.frame(0,0,0,0)
  for (i in 1:nrow(immu_pvalue_sig)) {#免疫细胞所在列
    c<-rcorr(immu1[,i],immu1[,mm],type ="spearman")
    cor1[i,1]<-c$r[2]
    cor1[i,2]<-c$P[2]
    cor1[i,3]<-colnames(immu1)[i]
    cor1[i,4]<-colnames(immu1)[mm]
  }
  cor<-rbind(cor,cor1)
}
colnames(cor)<-c("correlation","Pvalue","immu_cell","gene")
write.csv(cor, "output/7-immu_cor-cells-genes-Control.csv")
cor4<-cor#[cor$Pvalue<0.05,]
cor4<-na.omit(cor4)
cor4<-cor4[order(cor4$correlation,decreasing = F),]

CorPvalue <- cor4$Pvalue[which(cor4$Pvalue != 0)]
minCorPvalue <- min(CorPvalue)/2
for (i in c(1:length(cor4$Pvalue))){
  if (cor4$Pvalue[i] == 0){
    cor4$Pvalue[i]=minCorPvalue
  }
}

y=factor(cor4$gene)
x=factor(cor4$immu_cell)
a<-(-log10(cor4$Pvalue))
p<-ggplot(cor4,aes(x,y))+ geom_point(aes(size=a,color=correlation),shape=20)+
  scale_colour_gradient(high="#FF6A6A",low="#DADDFC")+
  scale_size_continuous(range = c(0,3))+
  # theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Low")+
  theme_classic()+labs(size="-log10Pvalue",x="",y="")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(colour = "black",size=6),
        legend.text = element_text(colour = "black",size=6),
        legend.title = element_text(colour = "black",size=6),
        legend.key.size = unit(0.2,"cm"),legend.spacing = unit(0.1, "cm"),legend.margin = ggplot2::margin(6,2,0,0),
        legend.background = element_blank(),plot.background = element_blank(),
        plot.margin = ggplot2::margin(6,0,0,0))+
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 2))
p
ggsave(plot = p,file="output/7-immu-cell-Control.pdf",width =8,height =6,units = "cm")
#####end#################################################################################
