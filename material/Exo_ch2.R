


#x = arima.sim(list(order = c(1,0,1),ar=0.9, ma = -0.8), n = 500,sd=0.1) to be explored
#########################
# Simulation of an arma(1,1)


set.seed=1007   # to keep
x = arima.sim(list(order = c(1,0,1),ar=5/6, ma = -1/6), n = 1000,sd=1)
xfit<-arima(x,order = c(1,0,1))
summary(xfit)

xhat<- fitted(xfit)
xres<- resid(xfit)
xfit$coef
acf2(xres)
Box.test(xres,lag=10,fitdf=2) # default=Box-Pierce
# Box.test(xres,lag=10,type="Box-Pierce",fitdf=2)
# Box.test(xres,lag=10,type="Ljung-Box",fitdf=2)


tsdiag(xfit)

# to see the fit quality on a part
xb<-as.ts(x[1:100])
xhatb<-as.ts(xhat[1:100])
ts.plot(cbind(xb,xhatb),type.plot="single",col=c(1,2))

# forecasting 20 ahead
x.pr = predict(xfit, n.headk=20)
x.pr
U = x.pr$pred + x.pr$se
L = x.pr$pred - x.pr$se
minx = min(x,L); maxx = max(x,U)
ts.plot(x, x.pr$pred, xlim=c(-100,550), ylim=c(minx,maxx)) 
lines(x.pr$pred, col="red", type="o")
lines(U, col="blue", lty="dashed")
lines(L, col="blue", lty="dashed")

model_auto<-auto.arima(x)
model_auto

#########################
# Model estimation for an ar(2) to model rec (Recruitment variable)
453 months over the years 1950-1987

# Using regression
acf2(rec, 48)     # will produce values and a graphic 
(regr = ar.ols(rec, order=2, demean=F, intercept=TRUE))  # regression
regr$asy.se.coef  # standard errors 
res0<- regr$resid
acf2(res0)
mean(fitted(regr)[-(1:2)])
mean(x=rec)

# Using Yule-Walker
rec.yw = ar.yw(rec, order=2)
rec.yw$x.mean  # = 62.26278 (mean estimate)
rec.yw$ar      # = 1.3315874, -.4445447  (parameter estimates)
sqrt(diag(rec.yw$asy.var.coef))  # = .04222637, .04222637  (standard errors)
rec.yw$var.pred  # = 94.79912 (error variance estimate)

(rec.pr$se)

rec.pr = predict(rec.yw, n.ahead=24)
U = rec.pr$pred + rec.pr$se
L = rec.pr$pred - rec.pr$se
minx = min(rec,L); maxx = max(rec,U)
ts.plot(rec, rec.pr$pred, xlim=c(1980,1990), ylim=c(minx,maxx)) 
lines(rec.pr$pred, col="red", type="o")
lines(U, col="blue", lty="dashed")
lines(L, col="blue", lty="dashed")
mean(rec)
rec.pr$pred[24]  # to be conmpared with the previous line result


# Using mle for ar
rec.mle = ar.mle(rec, order=2)
rec.mle$x.mean
rec.mle$ar
sqrt(diag(rec.mle$asy.var.coef))
rec.mle$var.pred

# perso: using general arima
rec_arima<- arima(rec,order=c(2,0,0))
summary(rec_arima)
rec_arima<- arima(rec,order=c(2,0,0),method="ML")
summary(rec_arima)
rec_arima<- arima(rec,order=c(2,0,0),method="CSS")
summary(rec_arima)
# the results are nearly the same. This comes from the fact
# in AR case things are simpler...

####################################################
####################################################
# Exo for dealing with seasonal arima
# load male_empl_W9 (employment figues for young men between 16 and
# 19 years of age, from jan.1971 to Dec 1981;  11 years, 132 data )
Man_empl<- scan("male_empl_W9.txt")  # when the file is in the w. dir.
Man_empl_ts<- ts (Man_empl,start=c(1971,1),frequency=12)
plot(Man_empl_ts)
plot(Man_empl,type="l")
acf2(Man_empl)

# plots for seasonality
monthplot(Man_empl_ts)
seasonplot(Man_empl_ts)

Man_diff1<- diff(Man_empl)
acf(Man_diff1,40)

Man_diff12<- diff(Man_empl,12)
acf(Man_diff12,40)

Man_diff<- diff(Man_diff12)
acf2(Man_diff,40)  # seems to be a seasonal ma(1)

Man_ma12<- arima(Man_empl,order=c(0,1,0),seasonal=list(order=c(0,1,1),period=12))
res_Man12<- Man_ma12$resid
acf2(res_Man12,40) # suggests also a ma(1) term


# trying model arima(0,1,1)(0,1,1)12
Man_ma112<-  arima(Man_empl,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12))
res_Man112<- Man_ma112$resid
acf2(res_Man112,40) #  seems OK



Man_fitted<- Man_empl-resid(Man_ma112)
plot.ts(cbind(Man_empl,Man_fitted),plot.type="single",col=c(1,2))




# prediction...something wrong?
Man_pred<-  predict(Man_ma112, n.ahead=12)
U = Man_pred$pred + Man_pred$se
L = Man_pred$pred - Man_pred$se
minx = min(Man_empl,L); maxx = max(Man_empl,U)
ts.plot(Man_empl[-(1:12)], Man_pred$pred, xlim=c(1971,1982), ylim=c(minx,maxx)) 
lines(Man_pred$pred, col="red", type="o")
lines(U, col="blue", lty="dashed")
lines(L, col="blue", lty="dashed")
mean(Man_empl)
Man_pred$pred[12]


##################################
##################################
# Trial for m.a.  via filter Spencer

Spenc<- c(-3,-6,-5,3,21,46,67,74,67,46,21,3,-5,-6,-3)
Spenc<- Spenc/320

tau_reduc <- t(Spenc) %*%(Spenc)
tau_reduc

x<- rep(1,100)
x
y<- filter(x,Spenc)
y


x<- seq(1:100)
x
y<- filter(x,Spenc)
y

x2<- x^2
x2
y<- filter(x2,Spenc)
y

x3<- x^3
x3
y<- filter(x3,Spenc)
y


t<- seq(1:100)

x<- cos(2*pi*t/4)
x
y<- filter(x,Spenc)
y

x<- cos(2*pi*t*2/4)
y<- filter(x,Spenc)
y
#############################################



####################################################
####################################################
# W4
# Monthly unemployed young women between ages 16 and 19 in the US
# from january 1961 to August 2002 (in thousands)
# 500 observations
library(forecast)
library(caschrono)
library(astsa)




Women_unpl<- scan(file = "D:/Data_R/Women_unpl.txt") # using absolute adress
Women_unpl<- scan("Women_unpl.txt")  # when the file is in the w. dir.
# here le W. directory is "D:/Data_R/"
Women_unpl<- scan(file.choose())  # using brooser


# saving and reading as a ts file
#write(Women_unpl, "Women_unpl_ts", ncolumns=1)
# xx <- read.ts("D:/Data_R/Women_unpl_ts",start=c(1961,1),frequency=12)
# xx <- Women_unpl<- scan(file = "D:/Data_R/Women_unpl_ts") # using absolute adress

par(mfrow=c(1,1))
plot(Women_unpl,type="l")
length(Women_unpl)



# convenient to have both versions std vector and series
# transform as a ts file and having dates in the x axis
Women_ts<- ts(Women_unpl,start=c(1961,1),frequency=12)
plot.ts(Women_ts)

# just to control seasonality. Nothing interesting
monthplot(Women_ts)
seasonplot(Women_ts)






# graphically fitting  a trend
plot(Women_unpl,type="l")
Wtime<- time(Women_unpl)
lines(lowess(Wtime,Women_unpl,0.05),col=2, lwd=2) # red
lines(lowess(Wtime,Women_unpl,0.10),col=3, lwd=2)  # green
lines(lowess(Wtime,Women_unpl,0.20),col=4, lwd=2) # blue
xfit<-lowess(Wtime,Women_unpl,0.20)[[1]]
fitmean<-lowess(Wtime,Women_unpl,0.20)[[2]]
plot(xfit,fitmean,type="l")
fitmean


# try a polynomial linear regression to fit the data (exercise)
# scrutenize the residuals
# could try also splines, kernel,...


# ACF and PACF
acf(Women_unpl,lag=40)    #  Stationarity???
pacf(Women_unpl)
xy.acfb(Women_unpl,lag=30)  # {caschrono}
acf2(Women_unpl)  # {astsa} seems compatibility problems




# we start fitting a stationary model arma(1,1)
arma11_mod<- arima(Women_unpl,order=c(1,0,1))
# the AR coefficient is  significative
# and very close to 1, which reinforce nonstationarity suspiscion!
sqrt(arma11_mod$var.coef[2,2]) # the MA coefficient is significative!
arma11_mod$coef/sqrt(diag(arma11_mod$var.coef))
# the AR coeff is very close to 1, which reinforce nonstationarity suspiscion!
# intercept significant
summary(arma11_mod)  #the informations already given above

# looking at residuals 
Women_res11<- (arma11_mod$resid)
plot(Women_res11)
acf(Women_res11,lag=10)
pacf(Women_res11,lag=10)
xy.acfb(Women_res11,lag=10)
Box.test(Women_res11,lag=10,fitdf=2) # default=Box-Pierce
Box.test(Women_res11,lag=10,type="Box-Pierce",fitdf=2)
Box.test(Women_res11,lag=10,type="Ljung-Box",fitdf=2)

# seems good but try this:
Box.test(Women_res11,lag=12,fitdf=2)
Box.test(Women_res11,lag=12,type="Box-Pierce",fitdf=2) 
Box.test(Women_res11,lag=12,type="Ljung-Box",fitdf=2)
acf(Women_res11,lag=12)
pacf(Women_res11,lag=12)
xy.acfb(Women_res11,lag=12)

# Do it also with lag=30
Box.test(Women_res11,lag=30,type="Box-Pierce",fitdf=2)
acf(Women_res11,lag=30)
pacf(Women_res11,lag=30)
xy.acfb(Women_res11,lag=30)
# So when taking not enough lags, we may get wrong findings

# differentiation
Women_diff <- diff(Women_unpl) # Y_t - Y_(t-1)
plot(Women_diff,type="l")
acf(Women_diff)   # reinforce nonstationarity suspicion
pacf(Women_diff)
acf2(Women_diff)
xy.acfb(Women_diff,lag=30,numer=FALSE,type="l")



# Test of nonstationarity : augmented unit root test of Dickey-Fuller
# from the library tseries
adf.test(Women_unpl,k=10) # p-value =0.8162 we don't reject "unit root"
adf.test(Women_unpl,k=40) # p-value =0.8162 we don't reject "unit root"

#### We go on  by fitting an order one integrated time series ARIMA(1,1,1) 
arima_model<- arima(Women_unpl,order=c(1,1,1))
arima_model$coef
arima_model$var.coef
arima_model$coef/sqrt(diag(arima_model$var.coef)) # the AR coefficient is non significative!
# the MA coefficient is  significative!
summary(arima_model)
arima_model$sigma2    # to be compared with the one for arma(1,1)...seems similar
Wom_res111<- ts(arima_model$resid)
var(Wom_res111)
plot.ts(Wom_res111)
Box.test(Wom_res111,lag=30,type="Box-Pierce",fitdf=2) # rejects the white noise hypothesis
acf(Wom_res111,lag=30) #exhibits saisonnality order 12
pacf(Wom_res111,lag=30) #exhibits saisonnality order 12
xy.acfb(Wom_res111,lag=30)  #  

Women_fit<- Women_unpl -Wom_res111
# Women_fit<- ts(arima_model$fitted)
# Women_fit<- predict(arima_model,3)
Women_unpl[1:50]
Women_fit[1:50]
plot.ts(cbind(Women_unpl[1:100],Women_fit[1:100]),plot.type="single",col=c(2,3))  
plot(Women_fit,type="l")

# difF2
Wom_diff_12 <- ts(diff(Women_unpl,lag=12)) # Y_t - Y_(t-12)
xy.acfb(Wom_diff_12,lag=60) # reinforce nonstationarity suspicion
# not so clear interpretation :seasonality on pacf

# diff1 + diff12  
Wom_diff1_12<- ts(diff(Women_diff,lag=12))
acf(Wom_diff1_12,lag=30)
pacf(Wom_diff1_12,lag=30)
xy.acfb(Wom_diff1_12,lag=50)


# We need to take seasonality into account
# coef AR1 put to 0
sarima_model<- arima(Women_unpl,order=c(0,1,1),
                     seasonal=list(order=c(1,0,1),period=12))
sarima_model
sWomen_res<- ts(sarima_model$res)
# plot(sWomen_res)
acf(sWomen_res,lag=30)
pacf(sWomen_res,lag=30)
sarima_model$coef
sarima_model$var.coef
sarima_model$coef/diag(sqrt(sarima_model$var.coef))
sarima_model$aic
sarima_model$bic
AIC(sarima_model)
AICc(sarima_model) # no doesn't work
BIC(sarima_model)
summary(sarima_model)

Box.test(sWomen_res,lag=30,type="Box-Pierce",fitdf=3)
Box.test(sWomen_res,lag=30,type="Ljung-Box",fitdf=3)
Box.test(sWomen_res,lag=30,type=c("Box-Pierce","Ljung-Box"),fitdf=3)

#automatic process best AIC but higher sigma^2
model_auto<-auto.arima(Women_unpl)
model_auto    # arima(0,1,1) !!!!


# Learning from a subset. We remove the last 12 observations
# which will be predicted with the model fitted on the 
# 488 first observations
length(Women_unpl)
Women_unpl2<-  ts(Women_unpl[1:488],start=c(1961,1),freq=12)

# sarima_model2<- arima(Women_unpl,order=c(0,1,1),
#                       seasonal=list(order=c(1,1,1),period=12))
# sWomen_res<- ts(sarima_model2$res)
# 
# sarima_model2<- arima(Women_unpl,order=c(0,1,1),
#                       seasonal=list(order=c(2,1,3),period=12))
# sWomen_res<- ts(sarima_model2$res)
# xy.acfb(sWomen_res,lag=50)
# Box.test(sWomen_res,lag=10,type="Box-Pierce",fitdf=6)


# model estimated on the basis of the 488 first observations (500-12)
# then the prevision of the last 12 observations is performed

# model selected by automatic procedure :model_0 
sarima_model2<- auto.arima(Women_unpl2)
summary(sarima_model2)
# ARIMA(0,1,1)(1,0,2)[12] selected
# ARIMA(2,1,2)(1,0,2)[12] with the complete dataset

sarima_model2<- arima(Women_unpl2,order=c(0,1,1),
                      seasonal=list(order=c(2,1,3),period=12))
sWomen_res<- ts(sarima_model2$res)
xy.acfb(sWomen_res,lag=50)
Box.test(sWomen_res,lag=10,type="Box-Pierce",fitdf=6)

sarima_model2<- arima(Women_unpl2,order=c(0,1,1),
                      seasonal=list(order=c(1,0,1),period=12))

sarima_model2<- arima(Women_unpl2,order=c(0,1,1))

sarima_model2<- arima(Women_unpl2,order=c(0,1,1),
                      seasonal=list(order=c(1,0,2),period=12))

sarima_model2<- arima(Women_unpl2,order=c(0,1,1),
                      seasonal=list(order=c(2,0,0),period=12))
sWomen_res<- ts(sarima_model2$res)
xy.acfb(sWomen_res,lag=50)
Box.test(sWomen_res,lag=20,type="Box-Pierce",fitdf=3)



sWomen_res2<- sarima_model2$res[14:488]
#Women_fit2<- fitted(sarima_model2)[14:488]
Women_fit2<- Women_unpl2[14:488] - sWomen_res2

Women_unpl22<- as.ts(Women_unpl2[14:488])
Women_unpl22<- ts(Women_unpl22,start=c(1962,2),frequency=12)
Women_end<- ts(Women_unpl[489:500],start=c(2001,9),freq=12)
Women_part2<- ts(Women_unpl[372:487],start=c(1992,1),freq=12)

num=length(Women_part2)
n.ahead<- 12
predd<-predict(sarima_model2,n.ahead)$pred
predd<- ts(predd,start=c(2001,9),freq=12)
predd


# y = ts(append(Women_unpl22, predd), start=c(1962,2), freq=12)
y = ts(append(Women_part2, predd), start=c(1992,1), freq=12)
# num=length(Women_unpl22)

y
Women_part2
predd


rmspe = sqrt(sarima_model2$sigma2)
#y[460:487]
upp = ts(y[(num+1):(num+n.ahead)]+2*rmspe, start=c(2001,9), freq=12)
low = ts(y[(num+1):(num+n.ahead)]-2*rmspe, start=c(2001,9), freq=12)
dev.new(height=2.25)
par(mar=c(3,3,1.5,1), mgp=c(1.6,.6,0),cex=.9)
plot.ts(y, main="", ylab="Women_unemployment", ylim=c(300,800), xlim=c(1992,2003))
lines(upp, lty=2);  lines(low, lty=2);  abline(v=2001.583, lty=3)
lines(Women_end)

predd2<- ts(predd,start=c(2001,9),frequency=12)
plot.ts(cbind(Women_end,predd2,upp,low),plot.type="single",ylim=c(250,1000))
xx = c(time(low), rev(time(upp)))
yy = c(low, rev(upp))
polygon(xx, yy, border=NA, col=gray(.5, alpha = .3))

Women_end
upp

fit1<- Women_fit2
fit2<- Women_fit2
fit3<- Women_fit2



length(sWomen_res2)
plot(sWomen_res2)
acf(sWomen_res2,lag=40)
pacf(sWomen_res2,lag=40)
sarima_model2$coef
sarima_model2$var.coef
sarima_model2$coef/sqrt(diag(sarima_model2$var.coef))

library(car)
t_stat(sarima_model2)

plot.ts(cbind(Women_unpl22,Women_fit2),plot.type="single",lty=1:2,type="l",col=c("black","red"))
plot.ts(cbind(Women_unpl22[100:200],Women_fit2[100:200]),plot.type="single",lty=1:2,type="l",col=c("black","red"))
err2<- sqrt(100*var(sWomen_res2))
err2


plot.ts(cbind(Women_unpl22[1:100],fit1[1:100],fit2[1:100],fit3[1:100]),plot.type="single",lty=1:4,type="l",col=c
        ("black","red","blue","green"))

### mod1 et mod2
plot.ts(cbind(Women_unpl22[100:200],fit1[100:200],fit2[100:200]),plot.type="single",lty=c(1,2,3),type="l",col=c
        ("black","red","blue"))

### mod1 et mod3
plot.ts(cbind(Women_unpl22[100:200],fit1[100:200],fit3[100:200]),plot.type="single",lty=c(1,2,4),type="l",col=c
        ("black","red","green"))

### mod2 et mod3
plot.ts(cbind(Women_unpl22[100:200],fit2[100:200],fit3[100:200]),plot.type="single",lty=c(1,3,4),type="l",col=c
        ("black","blue","green"))


library(forecast)   
(mod0<- auto.arima(Women_unpl))
mod0$var.coef
acf(mod0$res,lag=30)
pacf(mod0$res,lag=30)
err2<- sqrt(100*var(mod0$res))
err2

###########################################
###########################################





