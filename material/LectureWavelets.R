library(WaveletComp)
#The continuous wavelet transform
#Wavelet coefficients of a periodic signal
x = periodic.series(start.period = 50, length = 1000)
x = x + 0.2*rnorm(1000) # add some noise
#Computation of the wavelet transform
my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",loess.span = 0,dt = 1, dj = 1/250,lowerPeriod = 16,upperPeriod = 128,make.pval = TRUE, n.sim = 10)
#The spectral content of a signal is hidden in its wavelet coefficients
#periodic case
wt.image(my.w, color.key = "quantile", n.levels = 250,legend.params = list(lab = "wavelet power levels", mar = 4.7))
#A non periodic time series
t=1:1000
x = periodic.series(start.period = 20, end.period = 100, length = 1000)+0.2*rnorm(1000)
plot(t,x,'l')
#wavelet analyzis
my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 16,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))
wt.image(my.w, color.key = "quantile", n.levels = 250,legend.params = list(lab = "wavelet power levels", mar = 4.7))
#Reconstruction
reconstruct(my.w, plot.waves = FALSE, lwd = c(1,2),legend.coords = "bottomleft", ylim = c(-1.8, 1.8))
#A series with linear period
t=1:1000
x = periodic.series(start.period = 20, end.period = 100, length = 1000)+0.2*rnorm(1000)
plot(t,x,'l')
#wavelet analysis
my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 16,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))
#reconsruction
my.rec <- reconstruct(my.w)
x.rec <- my.rec$series$x.r # x: name of original series
#A series with two periods
t=1:1000
x1 <- periodic.series(start.period = 80, length = 1000)
x2 <- periodic.series(start.period = 30, length = 1000)
x <- x1 + x2 + 0.2*rnorm(1000)
plot(t,x,'l')
my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 16,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, n.levels = 250,legend.params = list(lab = "wavelet power levels") )
reconstruct(my.w,sel.period = 80, plot.waves = TRUE, lwd = c(1,2),legend.coords = "bottomleft")
my.w$Period[(my.w$Period > 79) & (my.w$Period < 81)]
#Average power
## Not run: 
## The following example is adopted from Liu et al., 2007:

series.length <- 6*128*24
x1 <- periodic.series(start.period = 1*24, length = series.length)
x2 <- periodic.series(start.period = 8*24, length = series.length)
x3 <- periodic.series(start.period = 32*24, length = series.length)
x4 <- periodic.series(start.period = 128*24, length = series.length)

x <- x1 + x2 + x3 + x4

plot(x, type = "l", xlab = "index", ylab = "", xaxs = "i",
     main = "hourly series with periods of 1, 8, 32, 128 days")

my.data <- data.frame(x = x)
my.wt <- analyze.wavelet(my.data, "x", loess.span = 0, 
                         dt = 1/24, dj = 1/20, 
                         lowerPeriod = 1/4,
                         make.pval = TRUE, n.sim = 10) 
wt.image(my.wt, color.key = "i", main = "wavelet power spectrum",
         legend.params = list(lab = "wavelet power levels (equidistant levels)"),
         periodlab = "period (days)", timelab = "time elapsed (days)")
# ,spec.time.axis = list(at = index.ticks, labels = index.labels))

## Plot of average wavelet power:
wt.avg(my.wt, siglvl = 0.05, sigcol = "red", periodlab = "period (days)")
#EMD for extracting signals
library(Rlibeemd)
x <- seq(0, 2*pi, length.out = 500)
signal <- sin(4*x)
intermittent <- 0.1 * sin(80 * x)
y <- signal * (1 + ifelse(signal > 0.7, intermittent, 0))
plot(x = x,y = y,type = "l")
# Decompose with EEMD
imfs <- eemd(y, num_siftings = 10, ensemble_size = 50, threads = 1)
plot(imfs)
#Separate high frequencies from low frequencies
HF_Comp<-rowSums(imfs[, 1:3])
ts.plot(HF_Comp)
LF_Comp<-rowSums(imfs[, 4:ncol(imfs)])
ts.plot(LF_Comp)

my.data <- data.frame(HF_Comp = HF_Comp)
my.w <- analyze.wavelet(my.data, "HF_Comp",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 16,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels for high frequencies") )

my.data <- data.frame(LF_Comp = LF_Comp)
my.w <- analyze.wavelet(my.data, "LF_Comp",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 16,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels for low frequencies") )

#Bivariate time series
sig1 <- periodic.series(start.period = 1*24, length = 24*96)
sig2 <- periodic.series(start.period = 2*24, length = 24*96)
sig3 <- periodic.series(start.period = 4*24, length = 24*96)
sig4 <- periodic.series(start.period = 8*24, length = 24*96)
sig5 <- periodic.series(start.period = 16*24, length = 24*96)
x1 <- sig1 + sig2 + 3*sig3 + sig4 + sig5 + 0.5*rnorm(24*96)
x2 <- sig1 + sig2 - 3*sig3 + sig4 + 3*sig5 + 0.5*rnorm(24*96)

my.data <- data.frame(x1 = x1, x2 = x2)
my.wc <- analyze.coherency(my.data, my.pair = c("x1","x2"),
                           loess.span = 0,
                           dt = 1/24, dj = 1/100,
                           lowerPeriod = 1/2,
                           make.pval = TRUE, n.sim = 10)

wc.image(my.wc, n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "", periodlab = "period (days)")

wc.avg(my.wc, siglvl = 0.01, sigcol = "red", sigpch = 20,
       periodlab = "period (days)")

wt.image(my.wc, my.series = "x1")
wt.image(my.wc, my.series = "x2")
