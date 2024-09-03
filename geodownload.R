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




###设置主题############################################################################################
library(ggplot2)
mytheme <- 
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"), # # 设置图片的边距，单位为厘米
        plot.title = element_text(size = 5,vjust = -2)) + #vjust = -2  图片title下移(hjust水平方向，vjust垂直方向)
  theme(panel.background = element_blank(), # 面板背景设为空白
        panel.grid = element_blank(), # 面板网格线设为空白
        panel.border = element_rect(fill = NA, size = 0.75 * 0.47)) + # 面板边框，填充透明，大小为 0.75*0.47
  theme(axis.line = element_line(size = 0.75 * 0.47), # 坐标轴线的粗细
        axis.text = element_text(size = 5, color = "black"),  # 坐标轴文本字体大小为 5，颜色为黑色
        axis.title = element_text(size = 5), # 坐标轴标题字体大小为 5
        axis.ticks = element_line(size = 0.75 * 0.47)) + # 坐标轴刻度线的粗细
  theme(legend.key = element_rect(fill = "white"),# 图例键的填充颜色为白色
        legend.key.size = unit(c(0.3, 0.3), "cm"), # 图例键的大小
        legend.title = element_text(size = 5),# 图例标题的字体大小
        legend.text = element_text(size = 5),# 图例文本的字体大小
        legend.margin = margin(0,0,-0.2,0,unit = "cm"),  # 图例的边距，向下移动 0.2 厘米
        legend.box.margin = margin(0,0,0,0,unit = "cm"),# 图例框的边距
        legend.box.spacing = unit(0, "cm"),#增加图注水平间距
        legend.background = element_blank(),  # 图例背景设为空白
        legend.spacing = unit(0, "cm"), # 图例之间的间距
        legend.box.background = element_blank())# 图例框的背景设为空白




group_col <- c("#D77071", "#6888F5")


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



dat_before_long <- combine_data %>% # 对合并前的数据进行长格式转换
  t() %>%  # 转置数据
  data.frame() %>% # 转换为数据框
  tibble::rownames_to_column("sample") %>%  # 将行名转换为 sample 列
  dplyr::mutate(group = gs) %>% # 添加 group 列，表示数据集来源
  dplyr::mutate(sample = factor(.$sample, levels = .$sample)) %>% # 将 sample 列转换为因子，并设置水平
  tidyr::gather(key = geneid, value, - c(sample, group)) # 将数据转换为长格式，geneid 为基因名称，value 为表达量

p_before <- 
  ggplot(dat_before_long, aes(sample, value)) + # 映射样本和表达量
  geom_boxplot(aes(fill = group), lwd = 0.1 *0.47, outlier.shape = NA) + # 添加箱线图，填充颜色根据 group 决定，不显示离群点
  scale_fill_manual(values = group_col) +# 手动设置填充颜色
  labs(title = "Before") + 
  mytheme + # 自己的主题
  theme(legend.position = 'top', # 图注位置（上）
        legend.direction = "horizontal", # 图注方向（水平）
        legend.title = element_blank(), # 隐藏图注标签
        axis.title.x = element_blank(), # 隐藏x轴标签
        axis.title.y = element_blank()) + # 隐藏y轴标签
  theme(
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) # 这句注释，下面两句运行，则不显示x轴文本（显示与否取决于样本数量，样本太多不建议展示）
    axis.text.x = element_blank(), # 隐藏x轴内容
    axis.ticks.x = element_blank() # 隐藏x轴刻度线
  ) + 
  coord_cartesian(clip = "off") # 解决上右框线看起来浅的问题

p_before ## 看下效果

ggsave("output/boxplot_before_combine.pdf", plot = p_before, units = "cm",width = 8,height = 4)


###合并绘制boxplot#############################################################################






dat_after_long <- adjusted_counts %>% 
  t() %>% 
  data.frame() %>% 
  tibble::rownames_to_column("sample") %>% 
  dplyr::mutate(group = gs) %>% 
  dplyr::mutate(sample = factor(.$sample, levels = .$sample)) %>% 
  tidyr::gather(key = geneid, value, - c(sample, group))

p_after <- 
  ggplot(dat_after_long, aes(sample, value)) + # 映射样本和表达量
  geom_boxplot(aes(fill = group), lwd = 0.1 *0.47, outlier.shape = NA) + # 以分组填充颜色，不要离群点, lwd柱子描边
  scale_fill_manual(values = group_col) +
  labs(title = "After") + 
  mytheme + # 自己的主题
  theme(legend.position = 'top', # 图注位置（上）
        legend.direction = "horizontal", # 图注方向（水平）
        legend.title = element_blank(), # 隐藏图注标签
        axis.title.x = element_blank(), # 隐藏x轴标签
        axis.title.y = element_blank()) + # 隐藏y轴标签
  theme(
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) # 这句注释，下面两句运行，则不显示x轴文本（显示与否取决于样本数量，样本太多不建议展示）
    axis.text.x = element_blank(), # 隐藏x轴内容
    axis.ticks.x = element_blank() # 隐藏x轴刻度线
  ) + 
  coord_cartesian(clip = "off") # 解决上右框线看起来浅的问题

p_after ## 看下效果

ggsave("output/boxplot_after_combine.pdf", plot = p_after, units = "cm",width = 8,height = 4)

##PCA#############################################################
### PCA矫正前----

count_matrix_pca <- prcomp(t(combine_data), scale. = T)
count_matrix_pcs <- data.frame(count_matrix_pca$x, group = gs)

## pca1,pca2的百分比
percentage <- round(count_matrix_pca$sdev / sum(count_matrix_pca$sdev) * 100, 2) # 计算每个主成分的贡献率
percentage <- paste(colnames(count_matrix_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))# 将主成分名称和贡献率拼接成字符串
group_col <- c("#D77071", "#6888F5")
pca_before <- ggplot(count_matrix_pcs, aes(x = PC1,y = PC2, color = group)) + # 创建 ggplot 对象，映射主成分和分组
  geom_point(aes(shape = group), size = 0.75) +# 添加散点图，点形状根据 group 决定
  stat_ellipse(level = 0.95, show.legend = F, geom = "polygon", aes(fill = group),  alpha = 0.2) +# 添加置信椭圆
  geom_hline(yintercept = 0, colour = "black", linetype = "longdash", size = 0.75 * 0.47) + # 添加水平虚线
  geom_vline(xintercept = 0, colour = "black", linetype = "longdash", size = 0.75 * 0.47) + # 添加垂直虚线
  scale_color_manual(values = group_col) +
  scale_fill_manual(values = group_col) +
  xlab(percentage[1]) +ylab(percentage[2]) + # 设置 x 轴和 y 轴标签
  labs(title = "Before") +
  mytheme +
  theme(panel.grid = element_line(colour = "grey90", size = 0.75 * 0.47, linetype = 1), # 设置网格线
        # plot.title = element_text(size = 7, hjust = 0.2, vjust = 0.5, margin = unit(c(0,0,0,0), "cm")),
        legend.position = 'top', # 图注位置（上）
        legend.direction = "horizontal",# 设置图例位置和方向
        legend.margin = margin(0, 0, -0.1, 0, "cm"),# 设置图例边距
        legend.title  = element_blank(), legend.key.size = unit(c(0.15, 0.15), "cm")) + # 隐藏图例标题，设置图例键大小
  coord_cartesian(clip = "off")

pca_before

ggsave("output/pca_before_combine.pdf", plot = pca_before, units = "cm",width = 4.5, height = 4.4) 

### PCA矫正后----

adjust_matrix_pca <- prcomp(t(adjusted_counts), scale. = T) ## scale. = T 归一化
adjust_matrix_pcs <- data.frame(adjust_matrix_pca$x, group = gs)

percentage <- round(adjust_matrix_pca$sdev / sum(adjust_matrix_pca$sdev) * 100, 2)
percentage <- paste(colnames(adjust_matrix_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))
group_col <- c("#D77071", "#6888F5")
pca_after <- ggplot(adjust_matrix_pcs, aes(x = PC1,y = PC2, color = group)) + 
  geom_point(aes(shape = group), size = 0.75) +
  stat_ellipse(level = 0.95, show.legend = F, geom = "polygon", aes(fill = group),  alpha = 0.2) +
  geom_hline(yintercept = 0, colour = "black", linetype = "longdash", size = 0.75 * 0.47) +
  geom_vline(xintercept = 0, colour = "black", linetype = "longdash", size = 0.75 * 0.47) +
  scale_color_manual(values = group_col) +
  scale_fill_manual(values = group_col) +
  xlab(percentage[1]) +ylab(percentage[2]) +
  labs(title = "After") +
  mytheme +
  theme(panel.grid = element_line(colour = "grey90", size = 0.75 * 0.47, linetype = 1), 
        # plot.title = element_text(size = 7, hjust = 0, vjust = 0.5, margin = unit(c(0,0,0,0), "cm")),
        legend.position = 'top', # 图注位置（上）
        legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.title  = element_blank(), legend.key.size = unit(c(0.15, 0.15), "cm")) + # 图注方向（水平）
  coord_cartesian(clip = "off")

pca_after

ggsave("output/pca_after_combine.pdf", plot = pca_after, units = "cm",width = 4.5, height = 4.4) 

