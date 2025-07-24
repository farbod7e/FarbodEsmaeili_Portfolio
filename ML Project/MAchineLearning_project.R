load("/home/sahooein/projects/def-ubcxzh/sahooein/finalproject.rda")
library(tree)
library(caret)
library(glmnet)
library(randomForest)
library(gbm)
library(xgboost)
library(dplyr)
library(ggplot2)
library(reshape)
library(ggpubr)

set.seed(1)
nrep=100

cv.cM = function(nfold = 5, model.data){
  fold.idx = sample(1:nfold, length(levels(factor(model.data$sample_id))), replace = T)
  cM = list(ClassificationTree = matrix(0,2,2))
  for(name in c("RandomForest","XGBoost","ElasticNet"))cM[[name]]=matrix(0,2,2)
  for(ii in 1:nfold){
    train.names<-levels(factor(model.data$sample_id))[fold.idx!=ii]
    train<-model.data[model.data$sample_id %in% train.names,]
    test<-model.data[which(model.data$sample_id %in% train.names == FALSE),]
    ### Classification tree
    test.detect<-model.data[which(model.data$sample_id %in% train.names == FALSE),]$detect_result
    tree.detection<-tree(detect_result ~ .,data = train[,-c(1:4)])
    tree.pred <- predict(tree.detection, test[,-c(1:4)],type = "class")
    # detection at virus level
    predictiondata<-data.frame(prediction = tree.pred, seg_id = test$seg_id, iso_id = test$iso_id, virus_name = test$virus_name, sample_id = test$sample_id)
    test.detection<-data.frame(prediction = test.detect, seg_id = test$seg_id, iso_id = test$iso_id, virus_name = test$virus_name, sample_id = test$sample_id)
    predictiondata<-predictiondata %>% 
      group_by(sample_id,virus_name) %>% 
      summarise(virus_name = dplyr::first(virus_name),
                seg_id = dplyr::first(seg_id),
                iso_id = dplyr::first(iso_id),
                prediction = any(as.logical(prediction)))
    test.detection<-test.detection %>% 
      group_by(sample_id,virus_name) %>% 
      summarise(virus_name = dplyr::first(virus_name),
                seg_id = dplyr::first(seg_id),
                iso_id = dplyr::first(iso_id),
                prediction = any(as.logical(prediction)))
    cM[[1]]<-cM[[1]]+table(pred = predictiondata$prediction, true_result = test.detection$prediction)
    
    ### Random forest
    rf.detection<- randomForest(detect_result~.,data=train[,-c(1:4)],mtry=10,ntree=50)
    yhat.detection<-predict(rf.detection,test[,-c(1:4,22)])
    predictiondata6<-data.frame(prediction = yhat.detection, seg_id = test$seg_id, iso_id = test$iso_id, virus_name = test$virus_name, sample_id = test$sample_id)
    predictiondata6<-predictiondata6 %>% 
      group_by(sample_id,virus_name) %>% 
      summarise(virus_name = dplyr::first(virus_name),
                seg_id = dplyr::first(seg_id),
                iso_id = dplyr::first(iso_id),
                prediction = any(as.logical(prediction)))
    cM[[2]]<-cM[[2]]+table(pred=predictiondata6$prediction, true_result=test.detection$prediction)
    
    ### XGBoost
    boost.detection<-gbm(as.numeric(detect_result)-1~., data = train[,-c(1:4)], distribution = "bernoulli", n.trees = 5000,interaction.depth = 4)
    
    yhat.detection<-predict(boost.detection,newdata = test[,-c(1:4)],n.trees=500,type="response")
    yhat.detection<-as.numeric(yhat.detection>0.6)
    predictiondata<-data.frame(prediction = yhat.detection, seg_id = test$seg_id, iso_id = test$iso_id, virus_name = test$virus_name, sample_id = test$sample_id)
    predictiondata<-predictiondata %>% 
      group_by(sample_id,virus_name) %>% 
      summarise(virus_name = dplyr::first(virus_name),
                seg_id = dplyr::first(seg_id),
                iso_id = dplyr::first(iso_id),
                prediction = any(as.logical(prediction)))
    cM[[3]]<-cM[[3]]+table(pred=predictiondata$prediction, true_result=test.detection$prediction)
    
    ### Elastic net
    cv_5 = trainControl(method = "cv", number = 5)
    det_elnet = train(detect_result ~ .,data = train[,-c(1:4)],method = "glmnet",trControl=cv_5)
    pred<-predict(det_elnet,test[,-c(1:4,22)])
    predictiondata<-data.frame(prediction = pred, seg_id = test$seg_id, iso_id = test$iso_id, virus_name = test$virus_name, sample_id = test$sample_id)
    predictiondata<-predictiondata %>% 
      group_by(sample_id,virus_name) %>% 
      summarise(virus_name = dplyr::first(virus_name),
                seg_id = dplyr::first(seg_id),
                iso_id = dplyr::first(iso_id),
                prediction = any(as.logical(prediction)))
    cM[[4]]<-cM[[4]]+table(pred=predictiondata$prediction, true_result=test.detection$prediction)
    
  }
  return(cM)
}
metriccal = function(m){
  f1score = (2*m[2,2])/(2*m[2,2]+m[1,2]+m[2,1])
  accuracy = (m[2,2]+m[1,1])/sum(m)
  recal = m[2,2]/(m[2,2]+m[2,1])
  precision = m[2,2]/(m[2,2]+m[1,2])
  return(list(f1score=f1score, accuracy=accuracy,recal=recal,precision=precision))
}

cMlist = list()
for(ii in 1:nrep){
  cM_eachrep = cv.cM(nfold = 5, model.data.newdata)
  cMlist[[ii]] = cM_eachrep
}

save(cMlist, file = "cMlist.Rdata")
load(file.choose())



acc= matrix(NA,nrep, 4)
dimnames(acc)=list(rep=c(1:nrep),model=c("ClassificationTree","RandomForest","XGBoost","ElasticNet"))
for(ii in 1:nrep){
  for(jj in 1:4){
    acc[ii,jj] = metriccal(cMlist[[ii]][[jj]])$accuracy
  }
}

recal= matrix(NA,nrep, 4)
dimnames(recal)=list(rep=c(1:nrep),model=c("ClassificationTree","RandomForest","XGBoost","ElasticNet"))
for(ii in 1:nrep){
  for(jj in 1:4){
    recal[ii,jj] = metriccal(cMlist[[ii]][[jj]])$recal
  }
}

precision= matrix(NA,nrep, 4)
dimnames(precision)=list(rep=c(1:nrep),model=c("ClassificationTree","RandomForest","XGBoost","ElasticNet"))
for(ii in 1:nrep){
  for(jj in 1:4){
    precision[ii,jj] = metriccal(cMlist[[ii]][[jj]])$precision
  }
}

f1= matrix(NA,nrep, 4)
dimnames(f1)=list(rep=c(1:nrep),model=c("ClassificationTree","RandomForest","XGBoost","ElasticNet"))
for(ii in 1:nrep){
  for(jj in 1:4){
    f1[ii,jj] = metriccal(cMlist[[ii]][[jj]])$f1score
  }
}

########Plots

par(mfrow=(c(2,2)))
par(mar=c(2,3,4,2))

m1 = melt(f1)
p1<-ggplot(m1, aes(x=model, y=value, color = model)) +
  geom_boxplot()+
  ggtitle("F1 Scores for Each Model")+
  xlab("Model")+ylab("Value")+
  theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust=1, size = 5),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 7))


m2 = melt(acc)
p2<-ggplot(m2, aes(x=model, y=value, color = model)) +
  geom_boxplot()+
  ggtitle("Accuracy for Each Model")+
  xlab("Model")+ylab("Value")+
  theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust=1, size = 5),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 7))

m3 = melt(recal)
p3<-ggplot(m3, aes(x=model, y=value, color = model)) +
  geom_boxplot()+
  ggtitle("Recall for Each Model")+
  xlab("Model")+ylab("Value")+
  theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust=1, size = 5),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 7))

m4 = melt(precision)
p4<-ggplot(m4, aes(x=model, y=value, color = model)) +
  geom_boxplot()+
  ggtitle("Precision for Each Model")+
  xlab("Model")+ylab("Value")+
  theme(axis.text.x = element_text(angle = 15, vjust = 0.5, hjust=1, size = 5),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 7))

ggarrange(p1, p2, p3, p4, common.legend = TRUE, legend = "bottom")

