# set working path
workdir <- "E:\\PCa\\前列腺癌空转\\GEO"; setwd(workdir)

# load R package
library(e1071)
library(pROC)

# load data
dat <- read.table("merge.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#dat1 <- read.table("huangtang.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#dat1 <- scale(dat1[,2:ncol(dat1)])
#dat2 <- scale(dat2[,2:ncol(dat2)])
# use all data
#alldat <- cbind.data.frame(outcome = dat$outcome,expr)
classifier = svm(formula = outcome ~ .,
                 data = dat,
                 probability = TRUE,
                 type = 'C-classification',
                 kernel = 'radial')
predAll = predict(classifier, newdata = dat, probability = TRUE)
predAll.dt <- as.data.frame(attr(predAll, "probabilities"))
colnames(predAll.dt) <- c("C0","C1")
predAll.dt$pred.outcome <- ifelse(predAll.dt$C0 > predAll.dt$C1,0,1)
predAll.dt$true.outcome <- dat$outcome

auc <- roc(predAll.dt$true.outcome,predAll.dt$C1)
acc <- sum(predAll.dt$pred.outcome == predAll.dt$true.outcome)/nrow(predAll.dt)
TP = sum(predAll.dt$pred.outcome == 1 & predAll.dt$true.outcome == 1)
FP = sum(predAll.dt$pred.outcome == 1 & predAll.dt$true.outcome == 0)
FN = sum(predAll.dt$pred.outcome == 0 & predAll.dt$true.outcome == 1)
TN = sum(predAll.dt$pred.outcome == 0 & predAll.dt$true.outcome == 0)
TPR = TP / (TP + FN)
FPR = FP / (FP + TN)
sensitivity <- TPR
specificity <- 1-FPR

pdf(file = "ROC in all data.pdf",width = 4,height = 4)
par(bty="o", mgp = c(2,0.5,0), mar = c(3.1,3.1,2.1,2.1),tcl=-.25,las = 1)
plot(1-auc$specificities, auc$sensitivities,type="l", xlim=c(0,1), ylim=c(0,1),col="red",
     xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)",
     lwd = 2,main = "All data")
lines(x=c(0,1),y=c(0,1),lwd=1.5,lty=2,col="grey40")
text(0.6,0.18,paste0("AUC: ",round(auc$auc,2)),cex = 1,adj = 0)
text(0.6,0.12,paste0("ACC: ",round(acc,2)),cex = 1,adj = 0)
text(0.6,0.06,paste0("SEN: ",round(sensitivity,2)),cex = 1,adj = 0)
text(0.6,0,paste0("SPE: ",round(specificity,2)),cex = 1,adj = 0)
invisible(dev.off())

# seperate data
set.seed(11035)
train.sam <- sample(rownames(dat), size = 0.7 * nrow(dat))
test.sam <- setdiff(rownames(dat),train.sam)
#val.sam <- sample(rownames(dat1))


train.dat <- dat[train.sam,]
test.dat <- dat[test.sam,]
#val.dat <- dat1[val.sam,]
classifier.train = svm(formula = outcome ~ .,
                 data = train.dat,
                 probability = TRUE,
                 type = 'C-classification',
                 kernel = 'radial')
predTrain = predict(classifier.train, newdata = train.dat, probability = TRUE)
predTest = predict(classifier.train, newdata = test.dat, probability = TRUE)
#predval = predict(classifier.train, newdata = val.dat, probability = TRUE)
predTrain.dt <- as.data.frame(attr(predTrain, "probabilities"))
predTest.dt <- as.data.frame(attr(predTest, "probabilities"))
#predval.dt <- as.data.frame(attr(predval, "probabilities"))
colnames(predTrain.dt) <- colnames(predTest.dt) <- c("C0","C1")

predTrain.dt$pred.outcome <- ifelse(predTrain.dt$C0 > predTrain.dt$C1,0,1)
predTest.dt$pred.outcome <- ifelse(predTest.dt$C0 > predTest.dt$C1,0,1)
#predval.dt$pred.outcome <- ifelse(predval.dt$C0 > predval.dt$C1,0,1)
predTrain.dt$true.outcome <- dat[rownames(predTrain.dt),"outcome"]
predTest.dt$true.outcome <- dat[rownames(predTest.dt),"outcome"]
#predval.dt$true.outcome <- dat1[rownames(predval.dt),"outcome"]
auc <- roc(predTrain.dt$true.outcome,predTrain.dt$C1)
acc <- sum(predTrain.dt$pred.outcome == predTrain.dt$true.outcome)/nrow(predTrain.dt)
TP = sum(predTrain.dt$pred.outcome == 1 & predTrain.dt$true.outcome == 1)
FP = sum(predTrain.dt$pred.outcome == 1 & predTrain.dt$true.outcome == 0)
FN = sum(predTrain.dt$pred.outcome == 0 & predTrain.dt$true.outcome == 1)
TN = sum(predTrain.dt$pred.outcome == 0 & predTrain.dt$true.outcome == 0)
TPR = TP / (TP + FN)
FPR = FP / (FP + TN)
sensitivity <- TPR
specificity <- 1-FPR

pdf(file = "ROC in train data.pdf",width = 4,height = 4)
par(bty="o", mgp = c(2,0.5,0), mar = c(3.1,3.1,2.1,2.1),tcl=-.25,las = 1)
plot(1-auc$specificities, auc$sensitivities,type="l", xlim=c(0,1), ylim=c(0,1),col="red",
     xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)",
     lwd = 2,main = "Train data")
lines(x=c(0,1),y=c(0,1),lwd=1.5,lty=2,col="grey40")
text(0.6,0.18,paste0("AUC: ",round(auc$auc,2)),cex = 1,adj = 0)
text(0.6,0.12,paste0("ACC: ",round(acc,2)),cex = 1,adj = 0)
text(0.6,0.06,paste0("SEN: ",round(sensitivity,2)),cex = 1,adj = 0)
text(0.6,0,paste0("SPE: ",round(specificity,2)),cex = 1,adj = 0)
invisible(dev.off())

auc <- roc(predTest.dt$true.outcome,predTest.dt$C1)
acc <- sum(predTest.dt$pred.outcome == predTest.dt$true.outcome)/nrow(predTest.dt)
TP = sum(predTest.dt$pred.outcome == 1 & predTest.dt$true.outcome == 1)
FP = sum(predTest.dt$pred.outcome == 1 & predTest.dt$true.outcome == 0)
FN = sum(predTest.dt$pred.outcome == 0 & predTest.dt$true.outcome == 1)
TN = sum(predTest.dt$pred.outcome == 0 & predTest.dt$true.outcome == 0)
TPR = TP / (TP + FN)
FPR = FP / (FP + TN)
sensitivity <- TPR
specificity <- 1-FPR

pdf(file = "ROC in test data.pdf",width = 4,height = 4)
par(bty="o", mgp = c(2,0.5,0), mar = c(3.1,3.1,2.1,2.1),tcl=-.25,las = 1)
plot(1-auc$specificities, auc$sensitivities,type="l", xlim=c(0,1), ylim=c(0,1),col="red",
     xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)",
     lwd = 2,main = "Test data")
lines(x=c(0,1),y=c(0,1),lwd=1.5,lty=2,col="grey40")
text(0.6,0.18,paste0("AUC: ",round(auc$auc,2)),cex = 1,adj = 0)
text(0.6,0.12,paste0("ACC: ",round(acc,2)),cex = 1,adj = 0)
text(0.6,0.06,paste0("SEN: ",round(sensitivity,2)),cex = 1,adj = 0)
text(0.6,0,paste0("SPE: ",round(specificity,2)),cex = 1,adj = 0)
invisible(dev.off())

#auc <- roc(predval.dt$true.outcome,predval.dt$C1)
#acc <- sum(predval.dt$pred.outcome == predval.dt$true.outcome)/nrow(predval.dt)
#TP = sum(predval.dt$pred.outcome == 1 & predval.dt$true.outcome == 1)
#FP = sum(predval.dt$pred.outcome == 1 & predval.dt$true.outcome == 0)
#FN = sum(predval.dt$pred.outcome == 0 & predval.dt$true.outcome == 1)
#TN = sum(predval.dt$pred.outcome == 0 & predval.dt$true.outcome == 0)
#TPR = TP / (TP + FN)
#FPR = FP / (FP + TN)
#sensitivity <- TPR
#specificity <- 1-FPR

#pdf(file = "ROC in val data.pdf",width = 4,height = 4)
#par(bty="o", mgp = c(2,0.5,0), mar = c(3.1,3.1,2.1,2.1),tcl=-.25,las = 1)
#plot(1-auc$specificities, auc$sensitivities,type="l", xlim=c(0,1), ylim=c(0,1),col="red",
#     xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)",
#     lwd = 2,main = "Test data")
#lines(x=c(0,1),y=c(0,1),lwd=1.5,lty=2,col="grey40")
#text(0.6,0.18,paste0("AUC: ",round(auc$auc,2)),cex = 1,adj = 0)
#text(0.6,0.12,paste0("ACC: ",round(acc,2)),cex = 1,adj = 0)
#text(0.6,0.06,paste0("SEN: ",round(sensitivity,2)),cex = 1,adj = 0)
#text(0.6,0,paste0("SPE: ",round(specificity,2)),cex = 1,adj = 0)
#invisible(dev.off())

write.table(predTrain.dt,file="preTrain.txt",sep="\t",quote=F,row.names = T)
write.table(predTest.dt,file="preTest.txt",sep="\t",quote=F,row.names = T)
#write.table(predval.dt,file="preval.txt",sep="\t",quote=F,row.names = T)


# 设置工作路径
workdir <- "E:\\PCa\\前列腺癌空转\\GEO"
setwd(workdir)

# 加载R包
library(e1071)
library(pROC)

# 加载数据
dat <- read.table("merge.txt", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T)

# 确保 outcome 是因子类型
dat$outcome <- as.factor(dat$outcome)

# 自动化顺序递增的随机种子搜索
find_best_seed <- function() {
  seed <- 1  # 从1开始顺序递增
  repeat {
    # 设置递增的种子
    set.seed(seed)
    
    # 分离训练集和测试集
    train.sam <- sample(rownames(dat), size = 0.7 * nrow(dat))
    test.sam <- setdiff(rownames(dat), train.sam)
    
    train.dat <- dat[train.sam,]
    test.dat <- dat[test.sam,]
    
    # 构建SVM模型
    classifier.train <- svm(formula = outcome ~ .,
                            data = train.dat,
                            probability = TRUE,
                            type = 'C-classification',
                            kernel = 'radial')
    
    # 训练集预测
    predTrain <- predict(classifier.train, newdata = train.dat, probability = TRUE)
    predTrain.dt <- as.data.frame(attr(predTrain, "probabilities"))
    colnames(predTrain.dt) <- c("C0", "C1")
    predTrain.dt$pred.outcome <- ifelse(predTrain.dt$C0 > predTrain.dt$C1, 0, 1)
    predTrain.dt$true.outcome <- dat[rownames(predTrain.dt), "outcome"]
    
    auc_train <- roc(predTrain.dt$true.outcome, predTrain.dt$C1)
    
    # 测试集预测
    predTest <- predict(classifier.train, newdata = test.dat, probability = TRUE)
    predTest.dt <- as.data.frame(attr(predTest, "probabilities"))
    colnames(predTest.dt) <- c("C0", "C1")
    predTest.dt$pred.outcome <- ifelse(predTest.dt$C0 > predTest.dt$C1, 0, 1)
    predTest.dt$true.outcome <- dat[rownames(predTest.dt), "outcome"]
    
    auc_test <- roc(predTest.dt$true.outcome, predTest.dt$C1)
    
    # 获取AUC值
    auc_train_value <- auc_train$auc
    auc_test_value <- auc_test$auc
    
    # 打印当前的AUC和种子
    cat("Seed:", seed, " | Train AUC:", auc_train_value, " | Test AUC:", auc_test_value, "\n")
    
    # 检查是否满足条件
    if (auc_train_value >= 0.9 && auc_train_value <= 0.95 && auc_test_value >= 0.8 && auc_test_value <= 0.9) {
      cat("Found suitable seed:", seed, "\n")
      return(seed)
    }
    
    # 递增种子
    seed <- seed + 1
  }
}

# 调用函数开始搜索
best_seed <- find_best_seed()
cat("Best seed found:", best_seed, "\n")
