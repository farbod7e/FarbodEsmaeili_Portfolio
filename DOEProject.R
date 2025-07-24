rm(list = ls())

set.seed(54)
sample(1:8,8)
set.seed(75)
sample(1:8,8)


#install.packages("conf.design")
require("conf.design")
A=rep(rep(c(-1,1),4),each=2)
B=rep(rep(c(rep(-1,2),rep(1,2)),2),each=2)
C=rep(rep(c(rep(-1,4),rep(1,4)),1),each=2)
X=matrix(c(rep(1,16),A,B,C),nrow=16,ncol=4)
Block = rep(c(1,2),8)
I=c(0.29,0.34,0.32,0.41,1.68,1.67,2.15,2.72,0.7,0.68,0.51,0.82,2.88,3.37,3.18,3.67)
V=c(84,84,81,84,79,80,83,81,84,82,82,85,79,83,77,83)/100
P=I*V
Data <- data.frame(P, A,B,C, Block)
res.lm<-lm(P~A*B*C+Block, data=Data)
summary(res.lm)
res.aov<-aov(P~A*B*C+Block, data=Data)
summary(res.aov)

#install.packages("daewr")
library(daewr)
fullnormal(coef(res.lm)[-1],alpha=.25)

#Projected model
res.lm.p<-lm(P~A*B*C+Block-A:B:C, data=Data)
summary(res.lm.p)
res.aov.p<-aov(P~A*B*C+Block-A:B:C, data=Data)
summary(res.aov.p)


#Final model - remove non-significant terms
res.lm.f<-lm(P~A*B*C+Block-A:B:C-A:C, data=Data)
summary(res.lm.f)
res.aov.f<-aov(P~A*B*C+Block-A:B:C-A:C, data=Data)
summary(res.aov.f)


#Residual Analysis
#Normality
residuals.f=res.aov.f$residuals
qqnorm(residuals.f, ylim=c(min(residuals.f),max(residuals.f)), main = "Normal Q-Q Plot for Residuals",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles- Modified",
       plot.it = TRUE, datax = FALSE)

qqline(residuals.f, datax = FALSE, distribution = qnorm)

#Test normality using Shapiro Wilks
shapiro.test(residuals.f)


#Check Variance
Fitted_values=res.aov.f$fitted.values
plot(Fitted_values,residuals.f,ylab="Residuals",xlab="Fitted Values", main = "Checking for Homogeneity of Variance")
abline(h=0, col = 2, lty = 2)

#Check Independence
plot(seq(1:length(residuals.f)),residuals.f,type="b",ylab="Residuals",xlab="Order", main = "Checking for Independence of Residuals")
abline(h=0, col = 2, lty = 2)

#Check Residuals vs Factors
par(mfrow=c(1,3))
plot(A,residuals.f,ylab="Residuals",xlab="Factor A: Potato Type", main = "Residuals VS. Factors", col = A+3)
abline(h=0, col = 2, lty = 2)

plot(B,residuals.f,ylab="Residuals",xlab="Factor B: Boiled", main = "Residuals VS. Factors", col = B+3)
abline(h=0, col = 2, lty = 2)

plot(C,residuals.f,ylab="Residuals",xlab="Factor C: State", main = "Residuals VS. Factors", col = C+3)
abline(h=0, col = 2, lty = 2)

par(mfrow=c(1,1))


###????
library(DescTools)
PostHocTest(res.aov.f<-aov(P~factor(A)*factor(B)*factor(C)+Block-A:B:C-A:C, data=Data)
, method= "hsd")

#plots

#install.packages("FrF2")

library(FrF2)
Data1 <- Data[Data$Block==1,1:4]
colnames(Data1) <- c("Power","PotatoType", "Boiled", "State")
P.lm.1 <- lm(Power ~ PotatoType*Boiled*State , data = Data1)
cubePlot(P.lm.1, "PotatoType", "Boiled", "State", size= 0.4,main = expression(paste(2^3, " Design of P(Watts) {Rep:1}")), round=4,y.margin.add=.1)
Data2 <- Data[Data$Block==2,1:4]
colnames(Data2) <- c("Power","PotatoType", "Boiled", "State")
P.lm.2 <- lm(Power ~ PotatoType*Boiled*State , data = Data2)
cubePlot(P.lm.2, "PotatoType", "Boiled", "State", size= 0.4, main = expression(paste(2^3, " Design of P(Watts) {Rep:2}")), round=4, y.margin.add=0.1
)


#Main effect plots
MEPlot(res.lm.f, main = c("Main effects plot for P(Watts)"))
par(mfrow=c(1,3))
interaction.plot(x.factor     = Data$A,
                 trace.factor = Data$B,
                 response     = Data$P,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o",
                 trace.label = "B: Boiled",
                 xlab = "A: Potato Type",
                 ylab = "Mean of P(Watts)",
                 )

interaction.plot(x.factor     = Data$A,
                 trace.factor = Data$C,
                 response     = Data$P,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o",
                 trace.label = "C: State",
                 xlab= "A: Potato Type",
                 ylab = "Mean of P(Watts)")

interaction.plot(x.factor     = Data$B,
                 trace.factor = Data$C,
                 response     = Data$P,
                 type="b",
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o",
                 trace.label = "C: State",
                 xlab= "B: Boiled",
                 ylab = "Mean of P(Watts)")
par(mfrow=c(1,1))
#install.packages("pid")
library(pid)

contourPlot(res.lm.f, 
            main="Contour plot", "A","B",
            colour.function=terrain.colors)
contourPlot(res.lm.f, 
            main="Contour plot", "A","C",
            colour.function=terrain.colors)
contourPlot(res.lm.f, 
            main="Contour plot", "B","C",
            colour.function=terrain.colors)

