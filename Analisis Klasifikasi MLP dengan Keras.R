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
