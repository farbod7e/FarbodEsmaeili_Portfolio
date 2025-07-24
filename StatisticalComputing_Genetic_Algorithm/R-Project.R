set.seed(0)
#Libraries
library(MTPS)
library(tidyverse)
library(hrbrthemes)
library(ggplot2)
library(viridis)
library(MASS)
library(tree)
library(glmnet)
library(caret)
library(rpart)
library(wesanderson)
library(reshape)



data(HIV)
str(XX)
str(YY)

layout(mat=matrix(1:2,2,1), heights=c(1,3))
par(mar=c(0.1,4.1,4.1,2.1))
image(1:nrow(YY), 1:ncol(YY),YY, xlab="sample index", ylab="drug index")
par(mar=c(4.1,4.1,0,2.1))
image(1:nrow(XX), 1:ncol(XX),XX, xlab="sample index", ylab="mutation index")

yBin <- as.matrix(YY)
cutoffs <- c(2,3,3,1.5,1.5) # cutoff value to be used to define drug resistance
for(ii in 1:5) yBin[,ii] <- (10^yBin[,ii] < cutoffs[ii])*1

layout(mat=matrix(1:2,2,1), heights=c(1,3))
par(mar=c(0.1,4.1,4.1,2.1))
image(1:nrow(yBin), 1:ncol(yBin),yBin, xlab="sample index", ylab="drug index")
par(mar=c(4.1,4.1,0,2.1))
image(1:nrow(XX), 1:ncol(XX),XX, xlab="sample index", ylab="mutation index")


#############################

dim(yBin)
dim(XX)

apply(yBin,2,table) #make it in table format for report 
#apply(XX,2,table)

XX.df = as.data.frame(XX)

# Cross Validation


geneticAlg = function(xx, response,
                      itr = 100, P = 20, m.rate = .01)
{
  n = length(response)
  m = length(xx[1,])
  r = matrix(0,P,1)
  phi = matrix(0,P,1)
  runs = matrix(0,P,m)
  runs.next = matrix(0,P,m)
  runs.aic = matrix(0,P,1)
  aics = matrix(0,P,itr)
  run = NULL
  best.aic = 0
  best.aic.gen = rep(0,itr)
  
  # INITIALIZES STARTING GENERATION, FITNESS VALUES
  set.seed(0) 
  for(i in 1:P){
    runs[i,] = rbinom(m,1,.5)
    run.vars = xx[,runs[i,]==1]
    g = glm(response~.,run.vars, family = binomial)
    runs.aic[i] = extractAIC(g)[2]
    aics[i,1] = runs.aic[i]
    if(i==1){
      run = runs[i,]
      best.aic = runs.aic[i]
    }
    else{
      if(runs.aic[i] < best.aic){
        run = runs[i,]
        best.aic = runs.aic[i]
      }
    }
  }
  r = rank(-runs.aic)
  phi = 2*r/(P*(P+1))
  best.aic.gen[1]=best.aic
  
  ## MAIN
  for(j in 1:itr-1){
    
    # BUILDS THE NEW GENERATION, SELECTING FIRST PARENT BASED ON
    # FITNESS AND THE SECOND PARENT AT RANDOM
    for(i in 1:P/2){
      parent.1 = runs[sample(1:P,1,prob=phi),]
      parent.2 = runs[sample(1:P,1),]
      pos = sample(1:(m-1),1)
      mutate = rbinom(m,1,m.rate)
      runs.next[i,] = c(parent.1[1:pos],parent.2[(pos+1):m])
      runs.next[i,] = (runs.next[i,]+mutate)%%2
      mutate = rbinom(m,1,m.rate)
      runs.next[P+1-i,] = c(parent.2[1:pos],parent.1[(pos+1):m])
      runs.next[P+1-i,] = (runs.next[P+1-i,]+mutate)%%2
    }
    runs = runs.next
    
    # UPDATES AIC VALUES, FITNESS VALUES FOR NEW GENERATION
    for(i in 1:P){
      run.vars = xx[,runs[i,]==1]
      g = glm(response~.,run.vars, family = binomial)
      runs.aic[i] = extractAIC(g)[2]
      aics[i,j+1] = runs.aic[i]
      if(runs.aic[i] < best.aic){
        run = runs[i,]
        best.aic = runs.aic[i]
      }
    }
    best.aic.gen[j+1]=best.aic
    r = rank(-runs.aic)
    phi = 2*r/(P*(P+1))
  }
  
  ## OUTPUT
  #run 		# BEST LIST OF PREDICTORS FOUND
  #best.aic 	# AIC VALUE
  return(list(best.run = run, best.aic = best.aic))
  ## PLOT OF AIC VALUES
  #plot(-aics,xlim=c(0,itr),ylim=c(50,425),type="n",ylab="Negative AIC",
  #	xlab="Generation",main="AIC Values For Genetic Algorithm")
  #for(i in 1:itr){points(rep(i,P),-aics[,i],pch=20)}
}
cv.cM = function(nfold, data.1, data.0, idx.1, idx.0)
{
  #copy code from above
  cM = list(GeneticAlg = matrix(0,2,2))
  for(name in c("Logistic","ElasticNet"))cM[[name]]=matrix(0,2,2)
  for(ii in 1:nfold){
    cat("fold",ii,"\n")
    data.train <<- rbind(data.1[idx.1!=ii,], data.0[idx.0!=ii,])
    data.test  <<- rbind(data.1[idx.1==ii,], data.0[idx.0==ii,])
    # Genetic Algorithm _ AIC
    start_time <- Sys.time()
    genetic.model = geneticAlg(xx = data.train[,-229], response = data.train$yy,
                  itr = 100, P = 20, m.rate = .01)
    data.train.genetic = cbind(data.train[,genetic.model$best.run==1], yy = data.train$yy)
    data.test.genetic = cbind(data.test[,genetic.model$best.run==1], yy = data.test$yy)
    genetic.fits = glm(yy~., data = data.train.genetic, family = binomial)
    end_time <- Sys.time()
    genetic.pred = predict(genetic.fits,data.test.genetic,type="response")
    genetic.pred = as.numeric(genetic.pred>.5)
    cM[["GeneticAlg"]] = cM[["GeneticAlg"]] + table(actual=data.test$yy,predicted=genetic.pred)
    f=function(timediff){runTime[[irt]][["GeneticAlg"]] <<- append(runTime[[irt]][["GeneticAlg"]] ,as.numeric(timediff))}
    f(end_time-start_time)
    #logistic
    start_time <- Sys.time()
    glm.fits=glm(yy~.,data = data.train ,family=binomial)
    end_time <- Sys.time()
    glm.probs=predict(glm.fits,data.test,type="response")
    glm.pred=rep(0,length(glm.probs))
    glm.pred[glm.probs>.5]=1
    cM[["Logistic"]] = cM[["Logistic"]] + table(actual=data.test$yy,predicted=glm.pred)
    f=function(timediff){runTime[[irt]][["Logistic"]] <<- append(runTime[[irt]][["Logistic"]] ,as.numeric(timediff))}
    f(end_time-start_time)
    #ElasticNet
    start_time <- Sys.time()
    cv_5 = trainControl(method = "cv", number = 5)
    myfit_elnet = train(
      yy ~ ., data = data.train,
      method = "glmnet",
      trControl = cv_5
    )
    myid <- order(myfit_elnet$results$Accuracy, decreasing = T)[1]
    mybest <- myfit_elnet$results[myid, ]
    elnet.mod <- glmnet(data.train[,-dim(data.train)[2]], data.train$yy, 
                        alpha = mybest$alpha, lambda = mybest$lambda, family = "binomial")
    end_time <- Sys.time()
    elnet.pred=predict(elnet.mod,as.matrix(data.test[,-dim(data.train)[2]]), s=mybest$lambda, type="class")
    cM[["ElasticNet"]] = cM[["ElasticNet"]] + table(actual=data.test$yy,predicted=elnet.pred)
    f=function(timediff){runTime[[irt]][["ElasticNet"]] <<- append(runTime[[irt]][["ElasticNet"]] ,as.numeric(timediff))}
    f(end_time-start_time)
  }
  #returning the list of cMs
  return(cM)
}  


nfold  = 5
cMlistDrugs= list(Drug1=NA,Drug2=NA,Drug3=NA,Drug4=NA,Drug5=NA)
runTime <<- list(Drug1=list(),Drug2=list(),Drug3=list(),Drug4=list(),Drug5=list())
f=function(){
  for(irt in 1:5){
    runTime[[irt]][["GeneticAlg"]] <<- vector()
    runTime[[irt]][["Logistic"]] <<- vector()
    runTime[[irt]][["ElasticNet"]] <<- vector()
  }
}
f()
for(eachDrug in 1:5){
  irt <<- eachDrug
  cMlist = list()
  yy = yBin[,eachDrug]
  
  datacom = data.frame(XX,yy=yy)
  datacom$yy=as.factor(datacom$yy)
  #split data into 0 and 1 for stratified fold sampling
  data.1 = datacom[datacom$yy==1,]
  data.0 = datacom[datacom$yy==0,]
  
  # Repeat CV for 50 times to investigate uncertainty of the comparison
  
  
  set.seed(0)
  nrep = 50
  idxmat.1 = replicate(nrep, sample(1:nfold, sum(datacom$yy==1), replace=T))
  idxmat.0 = replicate(nrep, sample(1:nfold, sum(datacom$yy==0), replace=T))
  
  
  for(jj in 1:nrep){cat(eachDrug,jj,"\n",sep=";");cMlist[[jj]]= cv.cM(nfold, data.1, data.0, idx.1=idxmat.1[,jj], idx.0=idxmat.0[,jj])}
  cMlistDrugs[[eachDrug]]=cMlist
}
save(cMlistDrugs, file = "/home/farbodes/projects/def-ubcxzh/farbodes/cMlistDrugs.RData")
save(runTime, file = "/home/farbodes/projects/def-ubcxzh/farbodes/runTime.RData")

load("cMlistDrugs.RData")
load("runTime.RData")

##criteria
nrep = 50
metriccal = function(m){
  f1score = (2*m[2,2])/(2*m[2,2]+m[1,2]+m[2,1])
  accuracy = (m[2,2]+m[1,1])/sum(m)
  recal = m[2,2]/(m[2,2]+m[2,1])
  precision = m[2,2]/(m[2,2]+m[1,2])
  return(list(f1score=f1score, accuracy=accuracy,recal=recal,precision=precision))
}
f1arr= array(NA, dim=c(nrep, 3, 5))
dimnames(f1arr)=list(rep=c(1:nrep),model=c("GeneticAlg", "Logistic", "ElasticNet"),drug=c(1:5))
for(drug in 1:5){
  for(ii in 1:nrep){
    f1arr[ii,1,drug] = metriccal(cMlistDrugs[[drug]][[ii]][["GeneticAlg"]])$f1score
    f1arr[ii,2,drug] = metriccal(cMlistDrugs[[drug]][[ii]][["Logistic"]])$f1score
    f1arr[ii,3,drug] = metriccal(cMlistDrugs[[drug]][[ii]][["ElasticNet"]])$f1score
  }
}

#f1arr
par(mfrow = c(1,5))
for(drug in 1:5){
  boxplot(f1arr[,,drug], border=1:5, main = paste0("Drug: ",drug,":", colnames(yBin)[drug]),
          ylab =ifelse(drug==1,"F1 Score",""),las=2, ylim=range(f1arr))
  if(drug==5)legend("topleft",paste0("Model ",1:3,": ", dimnames(f1arr)$model), bty="n", text.col=1:5, x.intersp=-.9)
}

#wilcxon test

for(drug in 1:5){
  for(model1 in 1:3){
    for(model2 in 1:3){
      if(model2>model1){
        wt=wilcox.test(f1arr[,c(model1, model2),drug]%*%c(1,-1))
        cat("Drug",drug,"model",colnames(f1arr)[model1],"model",colnames(f1arr)[model2],sep=":")
        print(wt)
      }
    }
  }
}





#############################################################################################
#accuracy
accarr= array(NA, dim=c(nrep, 3, 5))
dimnames(accarr)=list(rep=c(1:nrep),model=c("GeneticAlg", "Logistic", "ElasticNet"),drug=c(1:5))
for(drug in 1:5){
  for(ii in 1:nrep){
    accarr[ii,1,drug] = metriccal(cMlistDrugs[[drug]][[ii]][["GeneticAlg"]])$accuracy
    accarr[ii,2,drug] = metriccal(cMlistDrugs[[drug]][[ii]][["Logistic"]])$accuracy
    accarr[ii,3,drug] = metriccal(cMlistDrugs[[drug]][[ii]][["ElasticNet"]])$accuracy
  }
}
par(mfrow = c(1,5))
for(drug in 1:5){
  boxplot(accarr[,,drug], border=1:5, main = paste0("Drug: ",drug,":", colnames(yBin)[drug]),
          ylab =ifelse(drug==1,"Accuracy",""),las=2, ylim=range(accarr))
  if(drug==5)legend("topleft",paste0("Model ",1:3,": ", dimnames(accarr)$model), bty="n", text.col=1:5, x.intersp=-.9)
}

#############################################################################################
#recall
recarr= array(NA, dim=c(nrep, 3, 5))
dimnames(recarr)=list(rep=c(1:nrep),model=c("GeneticAlg", "Logistic", "ElasticNet"),drug=c(1:5))
for(drug in 1:5){
  for(ii in 1:nrep){
    recarr[ii,1,drug] = metriccal(cMlistDrugs[[drug]][[ii]][["GeneticAlg"]])$recal
    recarr[ii,2,drug] = metriccal(cMlistDrugs[[drug]][[ii]][["Logistic"]])$recal
    recarr[ii,3,drug] = metriccal(cMlistDrugs[[drug]][[ii]][["ElasticNet"]])$recal
  }
}
par(mfrow = c(1,5))
for(drug in 1:5){
  boxplot(recarr[,,drug], border=1:5, main = paste0("Drug: ",drug,":", colnames(yBin)[drug]),
          ylab =ifelse(drug==1,"Recall",""),las=2, ylim=range(recarr))
  if(drug==5)legend("topleft",paste0("Model ",1:3,": ", dimnames(recarr)$model), bty="n", text.col=1:5, x.intersp=-.9)
}
#############################################################################################
#precision
prearr= array(NA, dim=c(nrep, 3, 5))
dimnames(prearr)=list(rep=c(1:nrep),model=c("GeneticAlg", "Logistic", "ElasticNet"),drug=c(1:5))
for(drug in 1:5){
  for(ii in 1:nrep){
    prearr[ii,1,drug] = metriccal(cMlistDrugs[[drug]][[ii]][["GeneticAlg"]])$precision
    prearr[ii,2,drug] = metriccal(cMlistDrugs[[drug]][[ii]][["Logistic"]])$precision
    prearr[ii,3,drug] = metriccal(cMlistDrugs[[drug]][[ii]][["ElasticNet"]])$precision
  }
}
par(mfrow = c(1,5))
for(drug in 1:5){
  boxplot(prearr[,,drug], border=1:5, main = paste0("Drug: ",drug,":", colnames(yBin)[drug]),
          ylab =ifelse(drug==1,"Precision",""),las=2, ylim=range(prearr))
  if(drug==5)legend("topleft",paste0("Model ",1:3,": ", dimnames(prearr)$model), bty="n", text.col=1:5, x.intersp=-.9)
}


runTime_long = melt(runTime)
colnames(runTime_long) = c("runTime","Method","Drug")
head(runTime_long)
ggplot(runTime_long, aes(x = Drug, y = runTime, color = Method)) +  # ggplot function
  geom_boxplot()
avg_runtime=list()
for(drug in 1:5){
  avg_runtime[[drug]]=lapply(runTime[[drug]], mean)
}
avg_runtime_long = melt(avg_runtime)
colnames(avg_runtime_long) = c("runTime","Method","Drug")
ggplot(avg_runtime_long, aes(x = Drug, y = runTime, fill = Method)) +  # ggplot function
  geom_bar(stat = "identity",  position="dodge")+
  scale_fill_viridis(discrete=T, name="Model")+
  scale_fill_manual(values=wes_palette(n=3, name="Darjeeling2"))

