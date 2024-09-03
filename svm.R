rm(list = ls())
#机器学习##########################################################################
library(dplyr)
library(tibble)
library(readxl)

# 从文件中读取数据
Group <- data.table::fread("合并数据集分组-2.csv",check.names = F)
data_all <- data.table::fread("合并数据集表达谱-2.csv",check.names = F)%>% 
  tibble::column_to_rownames("V1")%>% as.data.frame()
gene_all <- read.csv("gene.csv")
gene_all <- gene_all$gene# 选择基因列



#3.SVM#############################################################################

set.seed(234)
source('msvmRFE.R')  
library(e1071)
library(parallel)
library(dplyr)

mpd <- Group
mdat <- data_all

#inter_s来源于分组比较图,只保留了具有差异的基因
inter_s <-gene_all

mdat <- mdat[inter_s,mpd$ID]
#mpd=group information #mdat=matrix #inter_s=25 genes matrix


#把mpd的分组列加入mdat矩阵中,3是分组文件中group所在列
input<-t(mdat) %>% as.data.frame()
input$ID<-rownames(input)
input <-inner_join(mpd,input,by="ID")
input<-as.data.frame(input)
rownames(input)<-input$ID
input<-input[,-1]
input$group<-as.factor(input$group)

##group 排序##
input<-input[order(input$group),]
dir.create("svm2")

svmRFE(input, k=10, halve.above=100)
nfold = 10
nrows = nrow(input)
str(input)
range(input[,2:14])
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100)
results
top.features = WriteFeatures(results, input, save=F)
head(top.features)
length(top.features$FeatureName)
#top.features有多少基因就是1：多少，例如有40个基因，最前面就是1：40，这一句基因越多，花费的时间越长。
featsweep = lapply(1:length(top.features$FeatureName), FeatSweep.wrap, results, input)

save(featsweep,file = "featsweep.RData")
# 画图
no.info = min(prop.table(table(input[,1])))

#提取错误率
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))


pdf("svm2/1-SVM_ErrorInters.pdf",width = 6,height = 5)
PlotErrors(errors, no.info=no.info) #查看错误率
dev.off()


pdf("svm2/2-SVM_AccuracyInters.pdf",width = 6,height = 5)
Plotaccuracy(1-errors,no.info=no.info) #查看准确率
dev.off()

which.min(errors)


write.table(top.features[1:which.min(errors), "FeatureName"],file="svm2/3-SVM-gene.txt",sep="\t",row.names = T,quote=F)
#####end#################################################################################

predict()