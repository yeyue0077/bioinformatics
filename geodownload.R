#1.数据合并#################################################################
rm(list = ls()) 
options(timeout = 1000000000)

## 创建文件夹
dir.create("raw") ## 存放下载的文件
dir.create("output") ## 存放生成的文件

library(GEOquery)
## 下载数据集压缩文件到raw文件夹中
gset2 <- getGEO('GSE46560', destdir = "raw", #在此处下载数据集矩阵文件
                AnnotGPL = F, ## 在此处不下载注释文件
                getGPL = F) ## 在此处不下载平台文件
gset3 <- getGEO('GSE48060', destdir = "raw", #在此处下载数据集矩阵文件
                AnnotGPL = F, ## 在此处不下载注释文件
                getGPL = F) ## 在此处不下载平台文件

## GSET1芯片平台整理----

library(magrittr)
library(limma)
library(dplyr)

## GSET2芯片平台整理----

dim(exprs(gset2[[1]])) 

### group分组文件----

pdata2 <- pData(gset2[[1]])
colnames(pdata2)
unique(pdata2$"source_name_ch1")
unique(pdata2$"platform_id")


group2 <- pdata2 %>%
  dplyr::select(ID="geo_accession", group0="source_name_ch1") %>%
  dplyr::mutate(group = dplyr::case_when( ## 针对分组信息，利用关键词定义分组
    grepl("with restenosis patient's Blood", group0) ~ "ISR", 
    grepl("without restenosis patient's Blood", group0) ~ "Control"
  )) %>% 
  dplyr::select(1,3) %>%
  dplyr::filter(group == "ISR" | group == "Control") %>%
  dplyr::arrange(desc(group))

data.table::fwrite(group2,"output/GSE46560_group.csv", row.names = F)

### 注释信息----
unique(pdata2$"platform_id")
gpl2 <- getGEO(unique(pdata2$"platform_id"), destdir = "raw", AnnotGPL = F)
gpl2 <- Table(gpl2)
colnames(gpl2)


gpl2 <- gpl2 %>% dplyr::select("ID", "Gene Symbol") %>%
  dplyr::filter(!grepl("///", .$"Gene Symbol")) %>%
  dplyr::filter(.$"Gene Symbol" != "") %>%
  dplyr::filter(!grepl("---", .$"Gene Symbol")) %>%
  dplyr::mutate(ID = as.character(ID)) %>%
  dplyr::rename(GENE_SYMBOL = "Gene Symbol")

### 探针替换----

data2 <- exprs(gset2[[1]])
ID <- rownames(data2) %>% as.character()

data2 <- data2[,group2$ID] %>%
  apply(2, as.numeric) %>%
  data.frame() %>%
  cbind(ID = ID, .) %>%
  dplyr::left_join(., gpl2, by = "ID") %>%
  aggregate(x = ., by = .$GENE_SYMBOL %>% list(), FUN = mean) %>%
  tibble::column_to_rownames("Group.1") %>%
  dplyr::select(-c("ID","GENE_SYMBOL"))

data.table::fwrite(data2,"output/GSE46560_Matrix.csv", row.names = T)

### 判断需不需要log2+1----

exprSet <- data2
range(exprSet) ## 查看最大最小值
exprSet = as.data.frame(exprSet)
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
## 如果LogC=T,就进行循环，如果存在小于0的数，就用非数填入
if (LogC) {
  for(i in 1:ncol(ex)){
    ex[which(ex[,i] < 0),i] <- NaN
  }
  exprSet <- log2(ex + 1) ## 将ex进行log2+1转化
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}
data2 <- exprSet
range(data2)

### 标准化----

boxplot(data2, outline = F, las = 2) ## 看看样本间是否有批次
data2_norm <- normalizeBetweenArrays(data2) %>% 
  data.frame()
data.table::fwrite(data2_norm,"output/GSE46560_Matrix_norm.csv", row.names = T)

##data3
dim(exprs(gset3[[1]]))
pdata3 <- pData(gset3[[1]]) ## 获取临床信息
colnames(pdata3) ## 方便后续复制列名
unique(pdata3$"characteristics_ch1.1") ## 打开pdata1，找出自己的分组信息列，并选择分组信息关键词
unique(pdata3$"platform_id") ## 芯片平台

group3 <- pdata3 %>% 
  dplyr::select(ID = "geo_accession", group0 = "characteristics_ch1.1") %>% ## geo_accession一般都是样本名，分组信息列就是刚才的
  dplyr::mutate(group = dplyr::case_when( ## 利用分组信息列中的关键词定义分组
    grepl("patient with recurrent events", group0) ~ "ISR",
    grepl("patient without recurrent events", group0) ~ "Control"
  )) %>% 
  dplyr::select(1,3) %>% 
  dplyr::filter(group == "ISR" | group == "Control") %>% ## 选择研究的样本
  dplyr::arrange(desc(group)) ## 对分组列排序，注意desc是降序，去掉就是升序，注意自己疾病和正常样本的顺序
data.table::fwrite(group3, "output/GSE48060_group.csv", row.names = F) ## 保存

### 注释信息----

unique(pdata3$"platform_id") ## 芯片平台
gpl3 <- getGEO(unique(pdata3$"platform_id"), destdir = "raw", AnnotGPL = F) ## 下载芯片平台
gpl3 <- Table(gpl3) ## 点开gpl3，找两列（探针id列，基因id列）
colnames(gpl3) ## 方便复制

gpl3 <- gpl3 %>% dplyr::select("ID", "Gene Symbol") %>% ## 选择探针名和基因名列
  dplyr::filter(!grepl("///", .$"Gene Symbol")) %>% ## 去除基因名中有///的
  dplyr::filter(.$"Gene Symbol" != "") %>% ## 去除基因名为空的
  dplyr::filter(!grepl("---", .$"Gene Symbol")) %>% ## 去除基因名中有---的
  dplyr::mutate(ID = as.character(ID)) %>% ## 探针列转字符
  dplyr::rename(GENE_SYMBOL = "Gene Symbol") ## 重命名基因列为GENE_SYMBOL

### 探针替换----

data3 <- exprs(gset3[[1]]) ## data1的表达矩阵
data3 <- data3[,group3$ID]
ID <- rownames(data3) %>% as.character() ## 提取矩阵探针列并转字符

data3 <- data3[,group3$ID] %>% ## 按group1排序样本
  apply(2, as.numeric) %>% ## 表达量数值化
  data.frame() %>% ## 转数据框
  cbind(ID = ID, .) %>% ## 补充探针列
  dplyr::left_join(., gpl3, by = "ID") %>%## 合并data1和gpl1
  aggregate(x = ., by = .$GENE_SYMBOL %>% list(), FUN = mean) %>% ## 相同ID取均值
  tibble::remove_rownames() %>% ## 移除行名
  tibble::column_to_rownames("Group.1") %>% ## 以Group.1为行名
  dplyr::select(-c("ID","GENE_SYMBOL")) ## 删除ID列和GENE_SYMBOL列

# data1 <-apply(data1,1,function(row){row-min(row)})
# data1 <-data1-min(data1)
### 判断需不需要log2+1----

exprSet <- data3
range(exprSet) ## 查看最大最小值
exprSet <- as.data.frame(exprSet)
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
## 如果LogC=T,就进行循环，如果存在小于0的数，就用非数填入
if (LogC) {
  for(i in 1:ncol(ex)){
    ex[which(ex[,i] < 0),i] <- NaN
  }
  exprSet <- log2(ex + 1) ## 将ex进行log2+1转化
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}
data3 <- exprSet
range(data3)
data.table::fwrite(data3, "output/GSE48060_Matrix.csv", row.names = T)
### 标准化----

boxplot(data3, outline = F, las = 2) ## 看看样本间是否有批次
data3_norm <- normalizeBetweenArrays(data3) %>% 
  data.frame()
data.table::fwrite(data3_norm, "output/GSE48060_Matrix_norm.csv", row.names = T)

save(data2,data3,group2,group3,file = "Data.Rdata")



##合并################################################################
library(sva)
library(ggplot2)
library("limma")
library('dplyr')
library(readxl)


B <- "GSE46560"
C <- "GSE48060"



gs <- factor( c(rep(B, length(group2$ID)),
                rep(C, length(group3$ID))), levels = c(B,C))# 创建一个因子变量，表示数据集来源

gene<-intersect(rownames(data2),rownames(data3))


data2<-data2[gene,]
data3<-data3[gene,]
combine_data<-cbind(data2,data3)
combine_group<-rbind(group2,group3)



batch <- c(rep("2", length(group2$ID)),
           rep("3", length(group3$ID)))
adjusted_counts<-normalizeBetweenArrays(combine_data)

adjusted_counts <- ComBat(adjusted_counts, batch = batch)

range(adjusted_counts)

adjusted_counts <- normalizeBetweenArrays(adjusted_counts)
boxplot(adjusted_counts, outline = F, las = 2)
range(adjusted_counts)

# save(data1,data2,data1_norm,data2_norm,group1,group2,combine_data,adjusted_counts,combine_group,file="数据下载矫正.RDATA")

write.csv(adjusted_counts,"output/合并数据集表达谱.csv")
write.csv(combine_group,"output/合并数据集分组.csv")
