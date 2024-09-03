library(readr)
library(dplyr) 
library("limma")
a <- read.csv("合并数据集表达谱.csv", row.names = 1)
b <- read.csv("合并数据集分组.csv",row.names = 1)

b <- b %>% arrange(desc(group))
group_order <- unique(b$ID) 
a_reordered <- a[, group_order]
write.csv(a_reordered,file="合并数据集表达谱-2.csv")
write.csv(b,file="合并数据集分组-2.csv")

##处理样本在前##保证上调下调正确##
modType=c(rep("ISR",10),rep("Control",32))

#创建线性模型。

design <- model.matrix(~0+factor(modType,levels = c("ISR","Control"),ordered = F))
colnames(design) <- c("ISR","Control")
#注意这里，在构建模型时，有一些数据集的处理组在前，对照组在后，且因为R里面默认的字符存储顺序为字母表顺序，
#如果处理组的首写字母的顺序在对照组的后面，构建模型时前三个又是处理组，
#那design出来的结果会上下调就完全相反。所以如果是这种情况的话，需要定义因子顺序。
#如：“modType=c(rep("DCM",40),rep("Control",40))
#design <- model.matrix(~0+factor(modType,levels = c("DCM","Control"),ordered = F)).”
rt <- a_reordered
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(ISR-Control,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',number=Inf)
write.table(allDiff, file = "deg_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#贝叶斯检验。
logFoldChange = 0.5
adjustP = 0.05

diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & P.Value < adjustP )), ]
#根据筛选条件筛选出差异基因。
diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & P.Value < adjustP )), ]
#筛选上调的基因
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & P.Value < adjustP )), ]
#筛选下调的基因
#如果想把这些基因的结果输出，可以用write.table函数write.table(diffUp,file="up.xls",sep="\t",quote=F,row.names=F)，
#要基因名字就把参数row.names设置为T.

write.table(diffSig,file="Sig.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.csv(diffUp,file="up.csv")
write.csv(diffDown,file="down.csv")
