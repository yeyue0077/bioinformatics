#
library(dplyr)  # 数据处理包
library(tibble)# 数据框操作包
library(readxl)# 读取 Excel 文件的包

# 从文件中读取数据
Group <- data.table::fread("合并数据集分组-2.csv",check.names = F)
data_all <- data.table::fread("合并数据集表达谱-2.csv",check.names = F)%>% 
  tibble::column_to_rownames("V1")%>% as.data.frame()
gene_all <- read.csv("gene.csv")
gene_all <- gene_all$gene# 选择基因列

# 数据处理
gene <- gene_all
#gene<-read.table("clipboard")
gene <- as.data.frame(gene)
exp <- data_all
exp <- as.data.frame(exp)
exp <- exp[gene$gene,]

# 转置数据
sle <- as.data.frame(t(exp))
expgroup<- Group
sle <- sle[expgroup$ID,]


# 添加分组信息
unique(expgroup$group)[1]
sle$Treat <- c(rep(1,nrow(expgroup[expgroup$group==unique(expgroup$group)[1],])),rep(0,nrow(expgroup[expgroup$group==unique(expgroup$group)[2],])))#手动添加分组信息，疾病组为1，对照组为0

library(glmnet)
library(pbapply)
iter.times <- 500 # 设置迭代次数，速度非常慢请耐心，例文是500次
# 
lasso_fea_list <- list()
lambda_all <- c('lambda.1se','lambda.min')
lambda_choose <- lambda_all[2]
list.of.seed <- 1:iter.times
lasso_fea_list <- pblapply(list.of.seed, function(x) {
  set.seed(list.of.seed[x])
  outcome <- sle$Treat
  xx <- as.matrix(sle[,gene_all])
  cvfit <- cv.glmnet(xx,outcome,family="binomial")
  fea <- rownames(coef(cvfit, s = lambda_choose))[coef(cvfit, s = lambda_choose)[, 1] != 0]
  if (is.element("(Intercept)", fea)) {
    lasso_fea <- sort(fea[-1])
  } else {
    lasso_fea <- sort(fea)
  }
  return(lasso_fea)
})

lasso_res <- NULL
for (i in 1:iter.times) {
  lasso_res <- rbind.data.frame(lasso_res,
                                data.frame(iteration = i,
                                           n.gene = length(lasso_fea_list[[i]]),
                                           genelist = paste0(lasso_fea_list[[i]], collapse = " | "),
                                           stringsAsFactors = FALSE),
                                stringsAsFactors = FALSE)
}

uniquelist <- unique(lasso_res$genelist)
uniquelab <- LETTERS[1:length(uniquelist)]
lasso_res$uniquelab <- NA
for (i in 1:length(uniquelist)) {
  lasso_res[which(lasso_res$genelist == uniquelist[i]),"uniquelab"] <- uniquelab[i]
}
lasso_res$label <- paste(lasso_res$n.gene,"genes",lasso_res$uniquelab,sep = "_") # 最终模型标签
# 
table(lasso_res$label)
tablelasso_res <- as.data.frame(table(lasso_res$label))

sel.lasso_res <- tablelasso_res[tablelasso_res$Freq==max(tablelasso_res$Freq),1]

droplevels(sel.lasso_res)
#选个数最多的那个
sel.iter <- lasso_res[which(lasso_res$label == droplevels(sel.lasso_res)),"iteration"][1]

#
set.seed(sel.iter) # 设置当前种子以复现该基因集,这里得到的是1，可以直接设置1来复现

outcome <- sle$Treat

xx <- as.matrix(sle[,gene_all])


dir.create("lasso")

# 绘制LASSO回归曲线图
lasso_fit=glmnet(xx,outcome,family="binomial",alpha=1,lambda = NULL)
pdf(file = "lasso/1-Likelihood.pdf",height = 5,width = 6)
plot(lasso_fit,xvar = "lambda",label = T)
dev.off()

#绘制LASSO回归10折交叉验证图
cvfit = cv.glmnet(xx,outcome,family="binomial")
pdf(file = "lasso/2-Lambda.pdf",height = 5,width = 6)
plot(cvfit)
dev.off()


#查看最佳lambda
cvfit$lambda.min
# 获取LASSO选出来的特征
myCoefs <- coef(cvfit, s="lambda.min")
fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[,1]!= 0]
if(is.element("(Intercept)", fea)) {
  lasso_fea <- fea[-1] # 去掉截距项并排序
  lasso_coef <- myCoefs@x[-1]; names(lasso_coef) <- lasso_fea
} else {
  lasso_fea <- fea
  lasso_coef <- myCoefs@x; names(lasso_coef) <- lasso_fea
}

write.table(lasso_fea,file="lasso/4-LASSO-gene.txt",sep="\t",row.names = T,quote=F)

####保存系数和变量###
lasso_coef_df <- data.frame(Gene=names(lasso_coef), Coefficient=lasso_coef) 
write.table(lasso_coef_df, file="lasso/4-LASSO-gene-coef.csv", row.names = F)
# 将数据框保存为TSV文件  
write.table(lasso_coef_df, file="lasso/4-LASSO-gene-coef.tsv", sep="\t", row.names=FALSE, quote=FALSE)

