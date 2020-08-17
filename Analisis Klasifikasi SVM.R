#### analisis klasifikasi menggunakan metode SVM ####
### 1. svm kernel linier 
# 1.1.tuning (best parameter)
set.seed(12345)
tuning_lin <- tune(svm, dataY~. , 
                   data = data_training, kernel="linear",
                   types = "C-clasification", ranges= list(
                     cost = c(0.1,0.01, 0.001,1 , 10 , 100)))
summary(tuning_lin)
# 1.2. model svm kernel linier
# model_lin <- svm(dataY~. , data_training, kernel = "linear")
model_lin <- svm(dataY~. , data_training, kernel = "linear", 
                 cost=0.01,scale=T,types = "C-clasification",
                 decision.value=T)
pred_train_lin <- predict(model_lin, data_training)
table(pred_train_lin, data_training$dataY)
pred_test_lin<-predict(model_lin, data_testting)
table(pred_test_lin, data_testting$dataY)
mean(pred_test_lin==data_testting$dataY)

### 2. svm kernel polynomial
# 2.1. tuning (best parameter)
set.seed(12345)
tuning_pol <- tune(svm, dataY~. , data = data_training, 
                   kernel="polynomial",types = "C-clasification",
                   ranges= list( cost = c(10,100,200,300,400,500)))
summary(tuning_pol)
# 2.2. model svm kernel polynomial
model_pol <- svm(dataY~., data_training, kernel="polynomial",
                 degree=3, cost=10, 
                 types="C-clasification", decision.value=T)
pred_train_pol <- predict(model_pol, data_training)
table(pred_train_lin, data_training$dataY)
pred_test_pol <- predict(model_pol, data_testting)
table(pred_test_pol, data_testting$dataY)
mean(pred_test_pol==data_testting$dataY)

### 3. svm kernel sigmoid
# 3.1 tuning (best parameter)
set.seed(12345)
tuning_sig <- tune(svm, dataY~. , data = data_training, 
                   kernel="sigmoid",
                   types = "C-clasification",
                   ranges= list(cost = c(
                     0.1,1,10,100,200,300), 
                     gamma=c(0.1,1,2,3,4,5)))
summary(tuning_sig)
# 3.2 model svm kernel sigmoid
model_sig <- svm(dataY~., data_training, 
                 kernel="sigmoid", cost=300,gamma=1,  
                 types="C-clasification", decision.value=T)
pred_train_sig <- predict(model_sig, data_training)
table(pred_train_sig, data_training$dataY)
pred_test_sig <- predict(model_sig, data_testting)
table(pred_test_sig, data_testting$dataY)
mean(pred_test_sig==data_testting$dataY)

### 4. svm kernel RBF
# 4.1. tuning (best parameter)
set.seed(12345)
tuning_RBF <- tune(svm, dataY~. , data = data_training, 
                   kernel="radial",types = "C-clasification",
                   ranges= list(
                     cost = c(0.1,0.01, 0.001,1 , 10 , 100),
                     gamma=c(1,2,3,4,5)))
summary(tuning_RBF)
# 4.2. model svm kernel RBF
model_rbf <- svm(dataY~., data_training, 
                 kernel="radial", cost=0.1, gamma=1,
                 types="C-clasification", decision.value=T)
pred_train_rad <- predict(model_rbf, data_training)
table(pred_train_rad, data_training$dataY)
pred_test_rad <- predict(model_rbf, data_testting)
table(pred_test_rad, data_testting$dataY)
mean(pred_test_rad==data_testting$dataY)
