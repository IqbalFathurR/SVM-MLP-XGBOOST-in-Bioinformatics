##### Analisis Klasifikasi menggunakan metode XGBoost ####
library(xgboost)
library(keras)
library(tensorflow)
library(caret)
library(plyr)
library(Matrix)
## convert data to matrix
dim(data_training)
trainx <- data_training[,1:80]
trainx <- as.matrix(trainx)
testx <-data_testting[,1:80]
testx<-as.matrix(testx)
target_trainy<-as.numeric(as.factor(data_training[,81]))
target_testy<-as.numeric(as.factor(data_testting[,81]))
target_trainy<-target_trainy-1
target_testy<-target_testy-1
data_train<-data.frame(trainx,target_trainy)
data_test<-data.frame(testx, target_testy)
dim(data_train)
dim(data_test)
trainm <- sparse.model.matrix(data_train$target_trainy~.-1,
                              data=data_train)
head(trainm)
train_label <-data_train[,"target_trainy"]
train_matrix <-xgb.DMatrix(data = as.matrix(trainm),
                           label=train_label)
testm<-sparse.model.matrix(data_test$target_testy~.-1, 
                           data = data_test)
test_label<-data_test[,"target_testy"]
test_matrix <-xgb.DMatrix(data = as.matrix(testm), 
                          label=test_label)
# tuning parameter
xgb_grid_par = expand.grid(nrounds = 500,
                           eta = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5),
                           max_depth = c(2,4,6,8),
                           gamma = c(1,2,4), 
                           subsample = c(0.25, 0.5, 0.75),
                           min_child_weight = c(1, 2, 3), 
                           colsample_bytree = 1)
xgb_control<- trainControl(number=5, verboseIter=TRUE)
library(data.table)
tr<-data.table(data_train, keep.rownames = F)
modFitxgb <- train(form=factor(target_trainy)~.,
                   data = tr,
                   method = "xgbTree",
                   metric = "Accuracy",
                   trControl = xgb_control,
                   tuneGrid = xgb_grid_par)
print(modFitxgb)
#Parameter  ie no of class
nc <- length(unique(train_label))
xgb_params <- list("objective"="multi:softprob",
                   "eval_metric"="mlogloss",
                   "num_class"=nc, method="xgbTree")
watchlist <- list(train=train_matrix,test=test_matrix)
#XGB Model
set.seed(1000)
bst_model <- xgb.train(params = xgb_params,
                       data = train_matrix,
                       nrounds = 500,
                       watchlist = watchlist,
                       eta=0.5,max_depth=2,
                       gamma=2,colsample_bytree=1,
                       min_child_weight=3, subsample=0.5)
#training & test error plot
e<-data.frame(bst_model$evaluation_log)
plot(e$iter, e$train_mlogloss, col='blue')
lines(e$iter, e$test_mlogloss, col='red')
min(e$test_mlogloss)
e[e$test_mlogloss == 0.692789,]
#feature importance
imp<-xgb.importance(colnames(train_matrix),
                    model = bst_model)
print(imp)
xgb.plot.importance(imp)
library(dplyr)
#prediction
p<-predict(bst_model, newdata = test_matrix)
predik<-matrix(p, nrow = nc, ncol = length(p)/nc) %>%
  t()%>%
  data.frame()%>%
  mutate(label=test_label, max_prob=max.col(., "last")-1)
table(Prediction= predik$max_prob, Actual=predik$label)
