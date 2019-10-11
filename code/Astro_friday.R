library("forecast")
library("caschrono")
library("tseries")
library("astsa")
library("ltsa")
library("expsmooth")
library("timeSeries")
library("fma")
library("WaveletComp")
library(TSA)
library(arfima)
library(fracdiff)
library(vars)

##########################################################
##########################################################
#######  long memory fractional differentiated time series
# example cours
lvarve<- log(varve)  # could use also Box-Cox transformation

# lambda <- BoxCox.lambda(varve)  
# lambda
# varve_bc <- BoxCox(varve,lambda)# data data varve transformed with the optimal lambda (0.397)

plot(lvarve,type="l")
acf2(lvarve)  # ACF appears slow decreasing
fr_lvarve<-arfima(lvarve)   
fr_lvarve
res<- residuals(fr_lvarve)[[1]]  
plot(res,type="l")
acf2(res)
# Fin example cours

# exercise cours
# Below just the beginning
dflvarve<-  diff(lvarve) 
plot(dflvarve,type="l")
acf2(dflvarve)   #
# could let think that arima(0,1,1) would be convenient
# adjust and see residuals


######################################
# ARMAX model:
# you can put a regresssion part in the model
# with exogeneous inputs

# regression with simple exogeneous inputs

temp_lh<- time(LakeHuron) - 1920
Lac_detr<- LakeHuron -fitted(lm(LakeHuron~temp_lh))
acf2(Lac_detr)

mod.lac<- arima(LakeHuron, order = c(2,0,0), xreg = temp_lh)

yhat<- fitted(mod.lac)
resid.lac<- residuals(mod.lac)
acf2(resid.lac)
plot.ts(cbind(LakeHuron,yhat),lty=1,las=1,
        plot.type="single", ylab="level",xlab="Year",col=c("black","red"))


# regression with simple exogeneous input and other time series
trend  = time(cmort) 
temp   = tempr - mean(tempr)
temp2  = temp^2
summary(fit <- lm(cmort~trend + temp + temp2 + part, na.action=NULL))

acf2(resid(fit), 52) # implies AR2
stats::acf2(resid(fit), 52) # if the previous doesn't work
plot(resid(fit),type="l")

fit_armax<-sarima(cmort, 2,0,0, xreg=cbind(trend,temp,temp2,part) )
res2<- resid(fit_armax$fit)
res2ts<- ts(res2)
acf2(res2ts,52)
plot(res2,type="l")
res2
Box.test(res2,lag=15,type="Box-Pierce",fitdf=2)

#######################################################

#######################################################
# vector autoregressive model

library(vars)
# First adjust a VAR(1) model with a linear trend (type ="both)
x = cbind(cmort, tempr, part)
summary(VAR(x, p=1, type="both"))  # "both" fits constant + trend


# compare several information criteria to choose the VAR order
# with still a linear trend
# Below we choose p=2 as suggested by BIC 
VARselect(x, lag.max=10, type="both")


summary(fit <- VAR(x, p=2, type="both"))
dev.new()
par(mfrow=c(3,3))
acf(resid(fit), 52)
stats::acf(resid(fit), 52)

serial.test(fit, lags.pt=12, type="PT.adjusted")

(fit.pr = predict(fit, n.ahead = 24, ci = 0.95))  # 4 weeks ahead
# dev.new()
fanchart(fit.pr)  # plot prediction + error

