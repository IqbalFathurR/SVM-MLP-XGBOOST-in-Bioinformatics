# Packages-packages
library(affy)
library(GEOquery)
library(Biobase)
library(simpleaffy)
library(affyPLM)
library(u133x3pcdf)
library(u133x3p.db)
library(genefilter)
library(AnnotationDbi)
## Input Data D:\Iqbal\Skripsi\Data\GSE18732_RAW
GSE18732<- list.celfiles("D:/Iqbal/Skripsi/Data/GSE18732_RAW",
                         full.names=T)
GSE18732.AFFY = ReadAffy(filenames=GSE18732)
GSE18732.AFFY
class(GSE18732.AFFY)
### get pheno data ###
GSET_GSE18732 <- getGEO(GEO="GSE18732",GSEMatrix =TRUE)
str(GSET_GSE18732)
data.GSE18732 <- exprs(GSET_GSE18732[[1]])
class(data.GSE18732)
dim(data.GSE18732)
data.GSE18732[1:5,1:4]
pheno_GSE18732 <- pData(phenoData(GSET_GSE18732[[1]]))
str(pheno_GSE18732)
varLabels(phenoData(GSET_GSE18732[[1]]))
dim(pheno_GSE18732)
View(pheno_GSE18732)
setwd('D://Iqbal/Skripsi/Draft Skripsi') 
#save pheno in csv file 
write.csv(pheno_GSE18732, file ="D://Iqbal/Skripsi/
          Draft Skripsi/pheno data/pheno_GSE18732.csv" )

#### analisis deskriptif ####
##BOXPLOT
library(dplyr)
dev.new(width=4+dim(GSET_GSE18732)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(GSET_GSE18732))
)/2),4,2,1))
title <- paste ("GSE18732", '/','GPL9486', 
                " selected samples", sep ='')
boxplot(data.gse_coba, boxwex=0.7, notch=T, main=title,
        outline=FALSE, las=2)
#Melihat Pheno
table(pheno_GSE18732$`gender:ch1`)
table(pheno_GSE18732$`age:ch1`)
#Barplot Gender
View(pheno_GSE18732)
des <-table(pheno_GSE18732$`gender:ch1`)
des <- as.data.frame(des) 
barplot(des$Freq, main='Gender',col=c('yellow','blue'),
        xlab= 'Gender', ylab='Pasien',
        names.arg= c("Female","Male"))
des<-count(pheno_GSE18732, "pheno_GSE18732$age")
colnames(pheno_GSE18732)
#Pie Chart Age
library(plyr)
des2<-count(pheno_GSE18732, "characteristics_ch1.11")
des2$characteristics_ch1.11 = as.character(gsub(
  "age:", "", des2$characteristics_ch1.11))
percent <- round(des2$freq/sum(des2$freq)*100)
des2 <- as.data.frame(des2)
lbls2 <- paste(des2$characteristics_ch1.11, percent,
               '%', sep=' ')
pie(des2$freq,main ='Age',label= lbls2, col=c(
  "red","orange","yellow","blue","green"),border="brown",
  cex=0.5)

#### pre processing data ####
library(affyPLM)
eset.dChip=threestep(GSE18732.AFFY,
                     background.method = "RMA.2", 
                     normalize.method="quantile",
                     summary.method="median.polish")
Ekspres.GSE18732 <- exprs(eset.dChip)
dim(eset.dChip)
#BOXPLOT
#set parameter and draw the plot
dev.new(width=4+dim(GSET_GSE18732)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(GSET_GSE18732))
)/2),4,2,1))
title <- paste ("GSE18732", '/','GPL9486', 
                " selected samples", sep ='')
boxplot(Ekspres, boxwex=0.7, notch=T, main=title, 
        outline=FALSE, las=2)

##### filtering #####
library(hgu133plus2.db)
filterdata <- nsFilter(eset.dChip, require.entrez =T,
                       var.func = IQR,
                       remove.dupEntrez = T,
                       var.cutoff = 0.5,
                       feature.exclude = "^AFFX")
log <-filterdata$filter.log
eset <- filterdata$eset
class(eset)
featureNames(eset) <- make.names(featureNames(eset))
#Melihat hasil filter (eset)
dim(eset) 
head(eset)
class(eset)
eset
# Membuat data microarray menjadi matrix/data frame
databaru <- exprs(eset)
dim(databaru)
class(databaru)
head(databaru)
View(databaru)
#Kelas Vektor ada 3(Normal, Impaired Glucose Tolerance dan Diabetes Melitus) 
datacl <- c(1,0,1,0,2,0,2,2,2,1,0,1,2,1,1,0,2,2,0,2,
            0,0,2,0,2,1,2,2,0,0,2,2,1,2,1,0,0,2,2,0,
            2,1,0,2,2,0,0,1,2,0,2,2,2,1,1,0,2,2,0,2,
            2,2,2,2,1,0,1,2,2,2,1,2,0,2,2,2,2,0,2,0,
            0,0,1,0,1,0,1,0,2,0,0,0,0,1,0,0,2,0,0,0,
            0,1,0,1,2,0,2,0,1,1,0,0,0,1,2,1,0,0)
length(datacl)
class(datacl)
datafac <- factor(datacl,levels=0:2, 
                  labels= c("NGT","IGT","DM"))
datafac
#### Filtering feature selection with multtest ####
library(multtest) 
datattest <- mt.teststat(databaru,datacl,test="f") 
class(datattest)
length(datattest)
qqnorm(datattest)
qqline(datattest)
length(datattest)
#Adjusting p-value (untuk melihat p-value yg sesuai dan tidak sesuai) 
rawp = 2 * (1 - pnorm(abs(datattest))) 
length(rawp)
prosedur = c("Bonferroni", "Holm", "Hochberg", "BH", "BY") 
adjusted = mt.rawp2adjp(rawp, prosedur) 
data <- adjusted$adjp[,] 
data1 <- data[order(adjusted$index), ]
head(data1)
dim(data1)
#mengambil kolom rawp
ffs <- data1[,1] 
class(ffs)
length(ffs)
ffs[1 : 10]
#Adjusting rawp 
datarawp <- data.frame(databaru, ffs) 
row.names(datarawp) <- row.names(databaru)
class(datarawp)
head(datarawp)
dim(datarawp) 
library(dplyr) 
datarawpfilterfinal <- subset(datarawp, ffs < 0.00000001)
rownames(datarawpfilterfinal)
class(datarawpfilterfinal)
dim(datarawpfilterfinal) 
head(datarawpfilterfinal)
## mendefinisikan data baru setelah filter
datadef <- datarawpfilterfinal[,1:118]
head(datadef)
dim(datadef)
colnames(datadef)
data_siap   = as.data.frame (t((datadef)))
head(data_siap)
dim(data_siap)
dataY = as.factor(datacl) 
dataY
data_use = as.data.frame(cbind(data_siap,dataY)) 
dim(data_use)
summary(data_use)
#Usepackages(analisis)
library(e1071)
library(pROC)
library(devtools)
library(caret)
NGT<-dplyr::filter(data_use, data_use$dataY==0)
IGT<-dplyr::filter(data_use, data_use$dataY==1)
DTM<-dplyr::filter(data_use, data_use$dataY==2)
#penentuan ukuran data traning dan testing
set.seed(123)
rasio=8/10
n<-round(nrow(NGT)*rasio)
sampel_NGT<-sample(1:nrow(NGT),n)
m<-round(nrow(IGT)*rasio)
sampel_IGT<-sample(1:nrow(IGT),m)
x<-round(nrow(DTM)*rasio)
sampel_DTM<-sample(1:nrow(DTM),x)
#data training dan testing
#training
training_NGT=NGT[sampel_NGT,]
dim(training_NGT)
training_IGT=IGT[sampel_IGT,]
dim(training_IGT)
training_DTM=DTM[sampel_DTM,]
dim(training_DTM)
#testing
testing_NGT=NGT[-sampel_NGT,]
dim(testing_NGT)
testing_IGT=IGT[-sampel_IGT,]
dim(testing_IGT)
testing_DTM=DTM[-sampel_DTM,]
dim(testing_DTM)
data_training=rbind(training_NGT,training_IGT,training_DTM)
data_testting=rbind(testing_NGT,testing_IGT,testing_DTM)
dim(data_training)
dim(data_testting)
View(data_training)
#save data training dan testing .csv
setwd("D://Iqbal/Skripsi/Draft Skripsi")
write.csv(data_training, 
          file ="D://Iqbal/Skripsi/Draft Skripsi/data skripsi/train80gen.csv" )
write.csv(data_testting, 
          file ="D://Iqbal/Skripsi/Draft Skripsi/data skripsi/test80gen.csv" )
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

##### Analisis Klasifikasi menggunakan metode MLP ####
library(keras)
library(tensorflow)
library(mlbench)
dim(data_training)
dim(data_testting)
#rubah data ke matrix
trainx <- data_training[,1:80]
trainx <- as.matrix(trainx)
testx <-data_testting[,1:80]
testx <- as.matrix(testx)

target_trainy<-as.numeric(as.factor(data_training[,81]))
target_testy<-as.numeric(as.factor(data_testting[,81]))
target_trainy<-target_trainy-1
target_testy<-target_testy-1
# one hot encoding
train_target <- to_categorical(target_trainy)
test_target <- to_categorical(target_testy)
# Epoch 200
# create sequential model
model.a<- keras_model_sequential()
model.a %>%
  layer_dense(units =20, activation = 'relu' ,
              input_shape = c(80)) %>%
  layer_dense(units =3) %>%
  layer_activation('softmax')
summary(model.a)
#compile
model.a %>% compile(loss= 'categorical_crossentropy',
                    optimizer='adam',
                    metrics=c('accuracy'))
#fit model
history.a<-model.a %>% fit(trainx, train_target, 
                           epochs= 200, batch_size=10, 
                           validation_split=0.2)
plot(history.a)
#evaluated model with data test
model1<- model.a %>%
  evaluate(testx, test_target)
#prediction & confussion matrix test data
prob<-model.a%>%
  predict_proba(testx)
pred.a<-model.a%>%
  predict_classes(testx)
tabel1<-table(Predicted=pred.a, Actual=target_testy)

## create sequential model
model.b<- keras_model_sequential()
model.b %>%
  layer_dense(units =40, activation = 'relu' ,
              input_shape = c(80)) %>%
  layer_dense(units =3) %>%
  layer_activation('softmax')
summary(model.b)
#compile
model.b %>% compile(loss= 'categorical_crossentropy',
                    optimizer='adam',
                    metrics=c('accuracy'))
#fit model
history.b<-model.b %>% fit(trainx, train_target, 
                           epochs= 200, batch_size=10,
                           validation_split=0.2)
plot(history.b)
#evaluated model with data test
model2<- model.b %>%
  evaluate(testx, test_target)
#prediction & confussion matrix test data
prob<-model.b%>%
  predict_proba(testx)
pred.b<-model.b%>%
  predict_classes(testx)
tabel2<-table(Predicted=pred.b, Actual=target_testy)

## create sequential model
model.c <- keras_model_sequential()
model.c %>%
  layer_dense(units =60, activation = 'relu',
              input_shape = c(80)) %>%
  layer_dense(units =3) %>%
  layer_activation('softmax')
summary(model.c)
#compile
model.c %>% compile(loss= 'categorical_crossentropy',
                    optimizer='adam',
                    metrics=c('accuracy'))
#fit model
history.c <-model.c %>% fit(trainx, train_target,
                            epochs= 200,batch_size=10, 
                            validation_split=0.2)
plot(history.c)
#evaluated model with data test
model3<- model.c %>%
  evaluate(testx, test_target)
#prediction & confussion matrix test data
prob<-model.c%>%
  predict_proba(testx)
pred.c<-model.c%>%
  predict_classes(testx)
tabel3<-table(Predicted=pred.c, Actual=target_testy)

## create sequential model
model.d<- keras_model_sequential()
model.d %>%
  layer_dense(units =80, activation = 'relu' ,
              input_shape = c(80)) %>%
  layer_dense(units =3) %>%
  layer_activation('softmax')
summary(model.d)
#compile
model.d %>% compile(loss= 'categorical_crossentropy',
                    optimizer='adam',
                    metrics=c('accuracy'))
#fit model
history.d<-model.d %>% fit(trainx, train_target, 
                           epochs= 200, batch_size=10,
                           validation_split=0.2)
plot(history.d)
#evaluated model with data test
model4<- model.d %>%
  evaluate(testx, test_target)
#prediction & confussion matrix test data
prob<-model.d%>%
  predict_proba(testx)
pred.d<-model.d%>%
  predict_classes(testx)
tabel4<-table(Predicted=pred.d, Actual=target_testy)

## create sequential model
model.e<- keras_model_sequential()
model.e %>%
  layer_dense(units =100, activation = 'relu',
              input_shape = c(80)) %>%
  layer_dense(units =3) %>%
  layer_activation('softmax')
summary(model.e)
#compile
model.e %>% compile(loss= 'categorical_crossentropy',
                    optimizer='adam',
                    metrics=c('accuracy'))
#fit model
history.e<-model.e %>% fit(trainx, train_target, 
                           epochs= 200, batch_size=10,
                           validation_split=0.2)
plot(history.e)
#evaluated model with data test
model5<- model.e %>%
  evaluate(testx, test_target)
#prediction & confussion matrix test data
prob<-model.e%>%
  predict_proba(testx)
pred.e<-model.e%>%
  predict_classes(testx)
tabel5<-table(Predicted=pred.e, Actual=target_testy)

## create sequential model
model.f<- keras_model_sequential()
model.f %>%
  layer_dense(units =100, activation = 'relu',
              input_shape = c(80)) %>%
  layer_dense(units =50, activation = 'relu') %>%
  layer_dense(units = 3) %>%
  layer_activation('softmax')
summary(model.f)
#compile
model.f %>% compile(loss= 'categorical_crossentropy',
                    optimizer='adam',
                    metrics=c('accuracy'))
#fit model
history.f<-model.f %>% fit(trainx, train_target,
                           epochs= 200,batch_size=10,
                           validation_split=0.2)
plot(history.f)

#evaluated model with data test
model6<- model.f %>%
  evaluate(testx, test_target)
#prediction & confussion matrix test data
prob<-model.f%>%
  predict_proba(testx)
pred.f<-model.f%>%
  predict_classes(testx)
tabel6<-table(Predicted=pred.f, Actual=target_testy)
## create sequential model
model.g<- keras_model_sequential()
model.g %>%
  layer_dense(units =80, activation = 'relu',
              input_shape = c(80)) %>%
  layer_dense(units =40) %>%
  layer_dense(units =3) %>%
  layer_activation('softmax')
summary(model.g)
#compile
model.g %>% compile(loss= 'categorical_crossentropy',
                    optimizer='adam',
                    metrics=c('accuracy'))
#fit model
history.g<-model.g %>% fit(trainx, train_target, 
                           epochs=200, batch_size=10,
                           validation_split=0.2)
plot(history.g)
#evaluated model with data test
model7<- model.g %>%
  evaluate(testx, test_target)
#prediction & confussion matrix test data
prob<-model.g%>%
  predict_proba(testx)
pred.g<-model.g%>%
  predict_classes(testx)
tabel7<-table(Predicted=pred.g, Actual=target_testy)
## create sequential model
model.h <- keras_model_sequential()
model.h %>%
  layer_dense(units =60, activation = 'relu',
              input_shape = c(80)) %>%
  layer_dense(units =30, activation = 'relu') %>%
  layer_dense(units =3) %>%
  layer_activation('softmax')
summary(model.h)
#compile
model.h %>% compile(loss= 'categorical_crossentropy',
                    optimizer='adam',
                    metrics=c('accuracy'))
#fit model
history.h <-model.h %>% fit(trainx, train_target, 
                            epochs=200,batch_size=10,
                            validation_split=0.2)
plot(history.h)
#evaluated model with data test
model8<- model.h %>%
  evaluate(testx, test_target)
#prediction & confussion matrix test data
prob<-model.h%>%
  predict_proba(testx)
pred.h<-model.h%>%
  predict_classes(testx)
tabel8<-table(Predicted=pred.h, Actual=target_testy)
## create sequential model
model.i<- keras_model_sequential()
model.i %>%
  layer_dense(units =40, activation = 'relu' ,
              input_shape = c(80)) %>%
  layer_dense(units =20, activation = 'relu') %>%
  layer_dense(units =3) %>%
  layer_activation('softmax')
summary(model.i)
#compile
model.i %>% compile(loss= 'categorical_crossentropy',
                    optimizer='adam',
                    metrics=c('accuracy'))
#fit model
history.i<-model.i %>% fit(trainx, train_target, 
                           epochs= 200, batch_size=10,
                           validation_split=0.2)
plot(history.i)
#evaluated model with data test
model9<- model.i %>%
  evaluate(testx, test_target)
#prediction & confussion matrix test data
prob<-model.i%>%
  predict_proba(testx)
pred.i<-model.i%>%
  predict_classes(testx)
tabel9<-table(Predicted=pred.i, Actual=target_testy)

######epoch 400
## create sequential model
model.e<- keras_model_sequential()
model.e %>%
  layer_dense(units =100, activation = 'relu' ,input_shape = c(80)) %>%
  layer_dense(units =3) %>%
  layer_activation('softmax')
summary(model.e)
#compile
model.e %>% compile(loss= 'categorical_crossentropy',
                    optimizer='adam',
                    metrics=c('accuracy'))
#fit model
history.e<-model.e %>% fit(trainx, train_target, epochs= 400,
                           batch_size=10, validation_split=0.2)
plot(history.e)
#evaluated model with data test
model5<- model.e %>%
  evaluate(testx, test_target)
#prediction & confussion matrix test data
prob<-model.e%>%
  predict_proba(testx)
pred.e<-model.e%>%
  predict_classes(testx)
tabel5<-table(Predicted=pred.e, Actual=target_testy)

## create sequential model
model.i<- keras_model_sequential()
model.i %>%
  layer_dense(units =40, activation = 'relu' ,input_shape = c(80)) %>%
  layer_dense(units =20, activation = 'relu') %>%
  layer_dense(units =3) %>%
  layer_activation('softmax')
summary(model.i)
#compile
model.i %>% compile(loss= 'categorical_crossentropy',
                    optimizer='adam',
                    metrics=c('accuracy'))
#fit model
history.i<-model.i %>% fit(trainx, train_target, epochs= 400,
                           batch_size=10, validation_split=0.2)
plot(history.i)
#evaluated model with data test
model9<- model.i %>%
  evaluate(testx, test_target)
#prediction & confussion matrix test data
prob<-model.i%>%
  predict_proba(testx)
pred.i<-model.i%>%
  predict_classes(testx)
tabel9<-table(Predicted=pred.i, Actual=target_testy)

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
