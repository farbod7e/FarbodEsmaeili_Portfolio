####read data####
library(readxl)
library(astsa)
Book1 <- read_excel("Book1.xlsx")
#View(Book1)
data = Book1$close

dates <- seq(as.Date("2018-02-01"), as.Date("2018-05-31"), by = "day")
t=ts(data, start = c(2018,2), frequency = 365)
time_series_data <- zoo::zoo(t, order.by = dates)

####plot data####
plot(time_series_data, ylab = "Daily Closing Price")


####model####
t = 1:length(data)
t2 = t ^ 2
model = lm(data ~ t + t2)
summary(model)
resids = residuals(model)
fits = fitted(model)

####plots####
par(mfrow=c(3,1))
hist(resids)
qqnorm(resids)
qqline(resids)
plot(resids ~ t, type = 'o')

####norm test####
shapiro.test(resids)
library(nortest)
ad.test(resids)

####runs test####
library(randtests)
runs.test(resids)

####variance test####
fac = factor(rep(1:10, each = 12))
library(car)
leveneTest(resids, group = fac)

####TS####
####variance test####
leveneTest(data, group = fac)

####box-cox####
boxCox(data ~ t, lambda = seq(from = -2, to = 2, by = .1))
data.new = log(data)
leveneTest(data.new, group = fac)
leveneTest(diff(data), group = fac[-1])
leveneTest(diff(data.new), group = fac[-1])

####mean test####
data.new.diff = diff(data.new)
fac.diff = fac[-1]
kruskal.test(data.new.diff, g = fac.diff)
leveneTest(data.new.diff, group = fac.diff)

####ACF-PACF####

par(mfrow=c(2,1))
acf(data.new.diff)
pacf(data.new.diff)
par(mfrow=c(1,1))
####model####
y = ts(data = data.new)
order = c(0, 1, 0)
library(forecast)
fit <- Arima(y = y, order = order)
plot(fit$fitted,col="red",lwd = 1, lty = 2, ylab= "Log of TS Data")
lines(y)
legend("topright", c("actual data", "fitted values"), col = c("black", "red"), lty = c(1,2))
fit <-
  auto.arima(
    y = data.new,
    d = 1,
    D = 0,
    start.p = 0,
    start.q = 0,
    max.p = 5,
    max.q = 5,
    seasonal = TRUE,
    start.P = 0,
    start.Q = 0,
    max.P = 3,
    max.Q = 3,
    ic = 'aicc',
    stepwise = FALSE,
    parallel = TRUE,
    num.cores = 6,
  )
fit
####forecast####
plot(forecast(fit, h = 10))
lines(fit$fitted,col="red",lwd = 1, lty= 2)
actual=c(55.83,56.83,56.43,53.61,54.46,53.76,53.51,52.29,51.97,44.39)
actual=log(actual)
lines(121:130,actual,col="black")
legend("topright", c("Point Estimates", "80% CI", "95% CI" , "actual data", "fitted values"), lty = c(0,0,0,1,2),
       col = c("blue", "grey", "grey",  "black", "red"),fill=c("blue","darkgrey","lightgrey",0,0),
       cex = 0.8)

data.frame(forecast(fit, h = 10),actual=actual)

sarima(data.new,0,1,0)
sarima(data.new,0,1,0)$ttable
acf2(fit$residuals)



