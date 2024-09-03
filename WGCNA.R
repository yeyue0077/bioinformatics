rm(list = ls())
##WGCNA###########################################################################
library("dplyr")
library("readxl")
library("limma")
library("WGCNA")


#读取文件,并对输入文件整理



rt_type<- data.table::fread("合并数据集分组-2.csv",header = T, check.names = F)
exp<-data.table::fread("合并数据集表达谱-2.csv",check.names = F) %>% as.data.frame()
rownames(exp)<-exp$V1
exp<-exp[,-1]
exp<-exp[,rt_type$ID]
group <- rt_type$group


ControlNum=length(group[group=="Control"])
ISRNum=length(group[group=="ISR"])

exp<-as.data.frame(exp)
range(exp)
#exp <- log2(exp+1)


# m.vars=apply(exp,1,var)#方差
# m.vars=apply(exp,1,mean)#均值
m.vars=apply(exp,1,sd)#计算每一行的标准差

#expro.upper=exp[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.1))[5]),]
#选择标准差大于其第五个十分位数的行（即选择变异度较大的行）
expro.upper = exp[order(apply(exp,1,sd), decreasing = T)[1:8000],]
gene<-data.table::fread("DEMRGs.csv")
gene<-gene$DEMRG
table(gene %in% rownames(expro.upper))#检查提取的基因名是否在之前筛选的高变异度基因中  


dim(expro.upper)
# 转置高变异度基因数据框，并准备进行后续分析 
datExpr=as.data.frame(t(expro.upper));

nGenes = ncol(datExpr)# 基因数量

nSamples = nrow(datExpr) # 样本数量

# 重新调整原始exp数据框的列顺序，以匹配datExpr,即提取之前筛选的高变异度基因矩阵
exp <- exp[colnames(datExpr),rt_type$ID]



# 构建一个矩阵，其中包含数值数据，并设置行名和列名  
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)# 对数据中的重复值进行平均处理（这是WGCNA包中的一个函数）
#data=log2(data+1)
###################################筛选基因方法
#data = data[order(apply(data,1,mad), decreasing = T)[1:8000],]
table(gene %in% rownames(data))
library(readxl)

#data=data[apply(data,1,sd)>0.01,]
###################################
datExpr0=t(data)


###检查缺失值
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}



################样本聚类########
sampleTree = hclust(dist(datExpr0), method ="average");
# dist()表示转为数值，method表示距离的计算方式，其他种类的详见百度
sizeGrWindow(12,9)
# 绘制样本树:打开一个尺寸为12 * 9英寸的图形输出窗口
# 可对窗口大小进行调整
# 如要保存可运行下面语句
pdf(file="1_sample_cluster.pdf",width=12,height=9);
par(cex = 0.6)	# 控制图片中文字和点的大小
par(mar =c(0,4,2,0))	# 设置图形的边界，下，左，上，右的页边距
plot(sampleTree, main ="Sample clustering to detectoutliers",sub="", xlab="", cex.lab = 1.5,
     cex.axis= 1.5, cex.main = 2)
# 参数依次表示：sampleTree聚类树，图名，副标题颜色，坐标轴标签颜色，坐标轴刻度文字颜色，标题颜色
# 其实只要包括sampleTree和图名即可
dev.off()


#########################删除离群样本##########
# 绘制阈值切割线
abline(h = 80,col="red"); # 高度15，颜色红色
# 确定阈值线下的集群
clust = cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
# 以高度15切割，要求留下的最少为10个
table(clust)	# 查看切割后形成的集合
# clust1包含想要留下的样本.
keepSamples = (clust==1)	# 将clust序号为1的放入keepSamples
datExpr = datExpr0[keepSamples, ]		
# 将树中内容放入矩阵datExpr中，因为树中剩余矩阵不能直接作为矩阵处理
nGenes =ncol(datExpr)	# ncol，crow分别表示提取矩阵的列数和行数
nSamples =nrow(datExpr)

###准备临床数据
traitData = data.frame(  
  Control = c(rep(0, ISRNum), rep(1, ControlNum)),  
  ISR = c(rep(1, ISRNum), rep(0, ControlNum))  
)
row.names(traitData)=colnames(data)
##去除离群的样本数据##
traitData <- traitData[rownames(traitData) !="GSM1167076",]
fpkmSamples = rownames(datExpr)
traitSamples =rownames(traitData)
sameSample=intersect(fpkmSamples,traitSamples)
datExpr=datExpr[sameSample,]
datTraits=traitData[sameSample,]

###分组数据也去除掉离群数据##
rownames(rt_type) <- rt_type$ID
rt_type <- rt_type[rownames(rt_type) !="GSM1167076",]
group <- rt_type$group
ControlNum=length(group[group=="Control"])
ISRNum=length(group[group=="ISR"])
###样品聚类
sampleTree2 = hclust(dist(datExpr), method = "average")

traitColors = data.frame(  
  V1 = c(rep("#FFFFFF", ISRNum), rep("#6888F5", ControlNum)), # ISR组颜色 (空白颜色)，控制组颜色 (低组)  
  V2 = c(rep("#D77071", ISRNum), rep("#FFFFFF", ControlNum))   # ISR组颜色 (高组)，控制组颜色 (空白颜色)  
)

pdf(file="2_sample_heatmap.pdf",width=20,height=8)
plotDendroAndColors(sampleTree2, traitColors, 
                    groupLabels = names(datTraits),
                    main = "ISR Sample dendrogram and trait heatmap")

dev.off()

##保存数据

nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits,file="step1_input.Rdata")

#########2挑选最佳阈值power#######
#rm(list = ls())  
#load("step1_input.Rdata")
##1.确定合适的软阈值：网络拓扑分析##
###power值散点图
enableWGCNAThreads()   #多线程工作
powers = c(1:20)       #幂指数范围1:20

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf(file="3_scale_independence-ISR.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");#颜色可以修改
abline(h=0.85,col="red") #颜色可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")#颜色可以修改
dev.off()

sft #查看最佳power值
softPower = sft$powerEstimate #最佳power值
#自己手动设置数量
# softPower =as.numeric(10)
save(sft, softPower, file = "step2_power_value.Rdata")



##################### 3.一步法构建加权共表达网络，识别基因模块 ####################

net <- blockwiseModules(datExpr,
                        power = softPower,
                        maxBlockSize = ncol(datExpr),
                        corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
                        networkType = "unsigned",
                        TOMType = "unsigned", 
                        minModuleSize = 40,    ##越大模块越少
                        mergeCutHeight = 0.25, ##越大模块越少
                        numericLabels = TRUE, 
                        saveTOMs = F,
                        verbose = 3
)
table(net$colors) 
# power: 上一步计算的软阈值
# maxBlockSize:计算机能处理的最大模块的基因数量(默认5000),16G内存可以处理2万个，
# 计算资源允许的情况下最好放在一个block里面。
# corType：计算相关性的方法；可选pearson(默认)，bicor。后者更能考虑离群点的影响。
# networkType:计算邻接矩阵时，是否考虑正负相关性；默认为"unsigned",可选"signed", "signed hybrid"
# TOMType：计算TOM矩阵时，是否考虑正负相关性；默认为"signed",可选"unsigned"。但是根据幂律转换的邻接矩阵(权重)的非负性，所以认为这里选择"signed"也没有太多的意义。
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，可存储起来供后续使用，
# mergeCutHeight: 合并模块的阈值，越大模块越少,一般为0.25
# minModuleSize: 每个模块里最少放多少个基因，设定越大模块越少
# 输出结果根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
# **0 (grey)**表示**未**分入任何模块的基因。


## 模块可视化，层级聚类树展示各个模块

# Convert labels to colors for plotting
moduleColors <- labels2colors(net$colors)
table(moduleColors)
# Plot the dendrogram and the module colors underneath
pdf("step3_genes-modules_ClusterDendrogram.pdf",width = 16,height = 12)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

save(net, moduleColors, file = "step3_genes_modules.Rdata")


#保存模块赋值和模块特征基因信息，以供后续分析。
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

#moduleLabels代表每个基因对应的模块序号（在基因名下面）
moduleLabels
#modulecolos代表每个基因对应的模块所属颜色（略去了基因名，但顺序同上图一样
moduleColors
#MEs是以基因模块为单位，各个样本在这个模块中的表达量（别管是怎么来的）
MEs



####################### 4.关联基因模块与表型 #####################################
rm(list = ls())  
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")



###模块与性状数据热图
# 明确基因和样本数
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# 用颜色标签重新计算MEs
# 按照模块计算每个module的MEs（也就是该模块的第一主成分）
# 按照下面的样式，得到的MEs0是将模块以颜色代表，各个样本中每个颜色的表达量，见下图
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
# 对给定的(特征)向量重新排序，使相似的向量(通过相关性测量)彼此相邻。
MEs = orderMEs(MEs0)	# 这时的MEs是将相似模块相邻后，即调整MEs0模块顺序后矩阵，这样在画模块性状关系图时，能够清晰捕捉特征
# 计算基因模块MEs 与 临床特征的相关性以及p值
# use 给出在缺少值时计算协方差的方法 
moduleTraitCor = cor(MEs, datTraits, use = "p")# 计算相关性系数
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)	# 计算P值
pdf(file="Module_trait-ISR.pdf",width=8,height=11)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
col1 <- colorRampPalette(colors = c("#6888F5","white","#D77071"),space="Lab")#定义颜色:低中高
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = col1(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,#相关性热图中的字体大小
               zlim = c(-1,1),
               main = paste("ISR Module-trait relationships"))
dev.off()


###输出每个模块的基因
dir.create("output")
for (mod in 1:nrow(table(moduleColors))){  
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0("output/9_",modules,"_genes.txt"),sep="\t",row.names=F,col.names=F,quote=F)
}
#检查模块p值
if(table(moduleTraitPvalue<0.05)["FALSE"]==2*nrow(moduleTraitPvalue)){
  print("No")
} else{
  print("Yes")
}
#提取p值小于0.05的模块
cccc<-as.data.frame(moduleTraitPvalue) 
cccc<-subset(moduleTraitPvalue,moduleTraitPvalue[,1]<0.05)
#创建gene_colour数据框,将基因探针和模块颜色合并成一个数据框。
gene_colour<-cbind(probes,moduleColors)
gene_colour<-as.data.frame(gene_colour)
#检查gene_colour中哪些探针（即基因）在gene向量中，并提取这些行。
gene_colour_tar<-gene_colour[gene_colour$probes %in% gene,]
#从cccc的行名中删除"ME"字符串。然后，从gene_colour_tar中筛选出模块颜色在cccc中的行。
cccc<-gsub("ME","",rownames(cccc))
gene_colour_tar1<-subset(gene_colour_tar,gene_colour_tar$moduleColors %in% cccc)


write.csv(gene_colour_tar1,"gene_colour_tar.csv")
# TCGA_tar<-exp[gene,]
# write.csv(TCGA_tar,"TCGA_tar.csv")


TOM = TOMsimilarityFromExpr ( datExpr, power = 14)
###输出每个模块的网络输入文件
for(mod in 1:nrow(table(moduleColors))){
  modules = names(table(moduleColors))[mod]
  probes = colnames(datExpr)
  inModule = (moduleColors == modules)
  modProbes = probes[inModule]
  modGenes = modProbes
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  outEdge = paste0("10_",modules , "_network.txt")
  outNode = paste0("10_",modules, "_nodes.txt")
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = outEdge,
                                 nodeFile = outNode,
                                 weighted = TRUE,
                                 threshold = 0.2,      ##############################最小距离设置
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
}
write.csv(moduleTraitCor,"相关系数.csv")
write.csv(moduleTraitPvalue,"p值.csv")



