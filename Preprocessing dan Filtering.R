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
