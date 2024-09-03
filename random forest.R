#4.随机森林##############################################################################
library(dplyr)  # 数据处理包
library(tibble)# 数据框操作包
library(readxl)# 读取 Excel 文件的包

# 从文件中读取数据
Group <- data.table::fread("合并数据集分组-2.csv",check.names = F)
data_all <- data.table::fread("合并数据集表达谱-2.csv",check.names = F)%>% 
  tibble::column_to_rownames("V1")%>% as.data.frame()
gene_all <- read.csv("gene.csv")
gene_all <- gene_all$gene# 选择基因列


#### Randomforest
# install.packages("randomForestSRC")
# install.packages("visNetwork")
# BiocManager::install("randomSurvivalForest")
library(randomForestSRC)
library(tidyverse) 
# BiocManager::install("htmltools")
library(survival)
library(randomForestSRC)
library(dplyr)



data <- data_all
# rownames(data)<-data[,1]
# data<-data[,-1]

# mgene<- readxl::read_xlsx("Input/gene.xlsx")
#mgene来源于分组比较图,只保留了具有差异的基因
mgene<-gene_all
bb<-intersect(rownames(data),mgene)
group <- Group
data1<-data[bb,group$ID]#只保留输入基因的表达矩阵，并且按分组表格中的顺序排序

write.csv(data1,"DEMRGs_expression.csv")
ferr_rt <- t(data1)#将表达矩阵转置，使得每一行代表一个样本
ferr_rt <- ferr_rt %>% 
  as.data.frame() %>% 
  mutate(group = group$group) %>% # 添加分组信息列
  select(group,everything()) # 将分组信息列移到第一列
write.csv(ferr_rt,"target_EXP.csv")
#判断是不是对照组在前
unique(ferr_rt$group)[1]
# 将对照组放在前面，其他组放在后面
ferr_rt<-rbind(ferr_rt[ferr_rt$group=="Control",],ferr_rt[ferr_rt$group!="Control",])
unique(ferr_rt$group)[1]# 再次查看分组的第一种类型，确保对照组在前
# 将对照组编码为1，其他组编码为0
ferr_rt$group <- ifelse(ferr_rt$group == unique(ferr_rt$group)[1], 1,0)
#unique(ferr_rt$group)[1]是参考组，对照组或者正常组

set.seed(234)
#install.packages ("randomForest")
library(randomForest) 
#ncol(ferr_rt)
ntree <- 1000 # 设置随机森林中决策树的数量
randomForest_output <- randomForest(x = ferr_rt[2:ncol(ferr_rt)],  # 特征矩阵，从第二列开始，因为第一列是分组信息  
                                    y = ferr_rt$group,   
                                    #type = "classification", # 明确指定进行分类  
                                    importance = TRUE,  # 计算特征重要性 
                                    ntree = ntree,  # 设置决策树数量 
                                    proximity = TRUE) # 计算样本之间的邻近度

dir.create("randomForest2")

pdf(file = "randomForest2/1-randomForest.pdf",width = 7,height = 6)
plot(randomForest_output, type = "l")
dev.off()
pdf(file = "randomForest2/2-randomForest-gene-output.pdf",width = 8,height = 6)
varImpPlot(randomForest_output, 
           n.var=length(mgene),#显示所有特征的重要性
           type = 2,# 显示 MeanDecreaseGini 重要性
           scale=FALSE, # 不进行标准化
           cex = 0.7)# 调整字体大小
dev.off()

rf_importances <- importance(randomForest_output, scale = FALSE) 

rf_importances2 <- rf_importances %>%
  as.data.frame() %>%# 转换为数据框
  arrange(desc(IncNodePurity))# 按IncNodePurity降序排序
IncNodePuritynumber <- 0.5# 设置IncNodePurity阈值
#筛选IncNodePurity节点纯度大于0.5的行（基因）
rf_importances3 <- rf_importances2[which(rf_importances2$IncNodePurity > IncNodePuritynumber),]


# 提取结果
rf_res <- ferr_rt %>%
  dplyr::select(1, one_of(rownames(rf_importances3)))
head(rf_res)
write.table(rf_res, file = paste0("randomForest2/3-randomforest_IncNodePurity-",IncNodePuritynumber,"-output.txt"), sep = "\t", col.names = NA,row.names = T, quote = F)
rf_res_gene <- as.data.frame(colnames(rf_res)[2:ncol(rf_res)])
colnames(rf_res_gene) <- "Random Forest"
write.csv(rf_res_gene, file = paste0("randomForest2/4-randomforest_gene-output.csv"), row.names = F, quote = F)
#随机森林部分结束，结束
#####end#################################################################################
