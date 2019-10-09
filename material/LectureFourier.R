library(astrodatR)
data("COUP_var")
plot(1:15145,COUP_var[,2],'l')
#Import the data and visualize
gx=scan("/home/marianne/Desktop/gx.dat")
t=1:5000
plot(t,gx[1:5000],pch=20)
#Classical time series
#simulate a sample of length 500 from a white noise.
x <- rnorm(500)
plot.ts(x)
#Simulate ARMA
x=arima.sim(model=list(ar=c(1,-0.25),ma=1),300,rand.gen=rnorm)
plot.ts(x)
acf(x)
#simulate farima
library(arfima)
x <- arfima.sim(1000, model = list(phi = .2, dfrac = .3))
plot.ts(x)
acf(x)
#Spectral density of classical time series 
#White noise
freq <- seq.int(0, 0.5, length.out = 500)
spec <- rep(1,500)
plot(freq, spec,'l')
#MA(1)
freq <- seq.int(0, 0.5, length.out = 500)
spec <- ((1 + 0.5*cos(2 * pi * freq ))^2 + (0.5*sin(2 * pi * freq ))^2)
plot(freq, spec,'l')
#AR(1)
freq <- seq.int(0, 1, length.out = 500)
spec <- 1/((1 + 0.9*cos(2 * pi * freq ))^2 + (0.9*sin(2 * pi * freq ))^2)
plot(freq, spec,'l')
#FARIMA(1,0.2,1)
freq <- seq.int(-0.7, 0.7, length.out = 500)
spec <- ((1-cos(2 * pi * freq))^2+(sin(2 * pi * freq))^2)^(-0.2)/((1 + 0.3*cos(2 * pi * freq ))^2 + (0.3*sin(2 * pi * freq ))^2)
plot(freq, spec,'l')
#Synthetic examples
t = 1:100
x1 = 2*cos(2*pi*t*6/100) + 3*sin(2*pi*t*6/100)
x2 = 4*cos(2*pi*t*10/100) + 5*sin(2*pi*t*10/100)
x3 = 6*cos(2*pi*t*40/100) + 7*sin(2*pi*t*40/100)
x = x1+x2+ x3 
par(mfrow=c(2,2))
plot.ts(x1, ylim=c(-10,10), main = expression(omega==6/100~~~A^2==13))
plot.ts(x2, ylim=c(-10,10), main = expression(omega==10/100~~~A^2==41))
plot.ts(x3, ylim=c(-10,10), main = expression(omega==40/100~~~A^2==85))
plot.ts(x, ylim=c(-16,16),main="sum")
#periodogram
P = abs(2*fft(x)/100)^2
f = 0:50/100
#Spectral ANOVA
#data
c1=cos(2*pi*t*6/100) 
s1=sin(2*pi*t*6/100)
c2=cos(2*pi*t*10/100) 
s2=sin(2*pi*t*10/100)
c3=cos(2*pi*t*40/100) 
s3=sin(2*pi*t*40/100)
#ANOVA
summary(lm(x~I(c1)+I(c2)+I(c3)+I(s1)+I(s2)+I(s3)))
plot(f, P[1:51], type="o", xlab="frequency", ylab="periodogram")
anova(lm(x~c1+c2+c3+s1+s2+s3))
#Practical data
freq = 0:32768/65536
I = (4/65536)*abs(fft(gx)/sqrt(65536))^2
plot(freq[2:60000],I[2:60000],type='l',xlab='Frequency')
#Averaged periodgram
fish=scan("/home/marianne/Desktop/Enseignement/2019-2020/Astrostats/Fourier/recruit.dat")
par(mfrow=c(2,1))
ish.per = spec.pgram(fish, taper=0, log="no")
abline(v=1/4, lty="dotted")
k = kernel("daniell", 4)
fish.ave = spec.pgram(fish, k, taper=0, log="no")
abline(v=c(.25,1,2,3), lty=2)
# Repeat above lines using rec in place of soi on line 3
fish.ave$bandwidth
# 0.0649519 = reported bandwidth
fish.ave$bandwidth*(1/12)*sqrt(12)
# 0.01875 = Bw
#Spectral Analysis Tools and tapering
par(mfrow=c(3,1))
fish.per = spec.pgram(fish, taper=0, log="no")
abline(v=1/4, lty="dotted")
k = kernel("daniell", 4)
fish.ave = spec.pgram(fish, k, taper=0, log="no")
abline(v=c(.25,1,2,3), lty=2)
fish.taper = spec.pgram(fish, k, taper=0.5, log="no")
abline(v=c(.25,1,2,3), lty=2)
#Multivariate examples
par(mfrow=c(3,1))
mfdeaths.spc <- spec.pgram(ts.union(mdeaths, fdeaths), spans = c(3,3))
plot(mfdeaths.spc, plot.type = "coherency")
plot(mfdeaths.spc, plot.type = "phase")
#Unequally spaced data
library(lomb)
# ibex contains an unevenly sampled time series
data(ibex)
plot(ibex$hours,ibex$temp)
lsp(ibex[2:3],)
lsp(ibex$temp,times=ibex$hours,type='period',ofac=5)