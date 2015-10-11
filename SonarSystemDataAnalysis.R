#################################################################
# read sonar system data (bearings and bearing rates) from file #
#################################################################
SSdat <- read.table("PV_1p5_03.txt", dec = ".")
SSdat <- read.table("PV_1p5_05.txt", dec = ".")
SSdat <- read.table("PV_1p5_1.txt", dec = ".")
SSdat <- read.table("PV_2p5_03.txt", dec = ".")
SSdat <- read.table("PV_2p5_05.txt", dec = ".")
SSdat <- read.table("PV_2p5_1.txt", dec = ".")
SSdat <- read.table("PV_0p03_03.txt", dec = ".")
SSdat <- read.table("PV_0p03_05.txt", dec = ".")
SSdat <- read.table("PV_0p03_1.txt", dec = ".")
SSdat <- read.table("PV_0p1_03.txt", dec = ".")
SSdat <- read.table("PV_0p1_05.txt", dec = ".")
SSdat <- read.table("PV_0p1_1.txt", dec = ".")
SSdat <- read.table("PV_0p25_03.txt", dec = ".")
SSdat <- read.table("PV_0p25_05.txt", dec = ".")
SSdat <- read.table("PV_0p25_1.txt", dec = ".")
SSdat <- read.table("PV_0p5_03.txt", dec = ".")
SSdat <- read.table("PV_0p5_05.txt", dec = ".")
SSdat <- read.table("PV_0p5_1.txt", dec = ".")
SSdat <- read.table("PV_0p75_03.txt", dec = ".")
SSdat <- read.table("PV_0p75_05.txt", dec = ".")
SSdat <- read.table("PV_0p75_1.txt", dec = ".")
SSdat <- read.table("PV_1_03.txt", dec = ".")
SSdat <- read.table("PV_1_05.txt", dec = ".")
SSdat <- read.table("PV_1_1.txt", dec = ".")
SSdat <- read.table("PV_2_03.txt", dec = ".")
SSdat <- read.table("PV_2_05.txt", dec = ".")
SSdat <- read.table("PV_2_1.txt", dec = ".")
# view the input
# edit(SSdat)
# first line contains pure signal (model bearing)
Bm <- SSdat[1,]
Bm <- unlist(Bm)
plot(Bm, type = "l")
# model bearing rate
dBm <- diff(Bm)
dBm.mean <- mean(dBm); dBm.mean
plot(dBm, type = "l")
dBm
##############################################################
# lines 2-31 contain estimated bearing (sonar system output) #
##############################################################
Be <- SSdat[4,]
Be <- unlist(Be)
plot(Be, type = "l")
# seasonality analysis
Be.ts <- ts(Be, frequency = 10); Be.ts
Be.stl <- stl(Be.ts, "periodic"); Be.stl
Be.stl.seas <- Be.stl$time.series[,1]  # seasonality
Be.stl.trend <- Be.stl$time.series[,2] # trend
Be.stl.rem <- Be.stl$time.series[,3]   # remainder
plot(Be.stl.seas, type = "l")
plot(Be.stl.trend, type = "l")
plot(Be.stl.rem, type = "l")
# get noise signal
ksi <- Bm - Be
# analyse noise signal
plot(ksi, type = "l")
# mean
ksi.mean <- mean(ksi); ksi.mean
# mean root-square error
ksi.mrse <- sqrt(var(ksi)); ksi.mrse
# compute autocovariance
ksi.acf <- acf(ksi, lag.max = 3,type = "correlation")
ksi.acf <- acf(ksi, type = "correlation")
ksi.acf <- acf(ksi, type = "covariance")
ksi.acf <- acf(ksi, type = "partial")
ksi.acf
# spectral density
ksi.spec <- spectrum(ksi)
# cumulative periodogram
ksi.cpgram <- cpgram(ksi)
# k-step differences (lag = k)
zeta <- diff(ksi, lag = 1, differences = 1)
plot(zeta, type = "l")
zeta.acf <- acf(zeta, type = "correlation")
zeta.cpgram <- cpgram(zeta)
# ARIMA model for residuals
ksi.arima <- arima(ksi, order = c(10,0,0)); ksi.arima
cpgram(ksi.arima$resid, main = "ARIMA fit to residuals")
# ARIMA model with Akaike's information criterion (AIC) model selection
ksi.ar <- ar(ksi, order.max = 10); ksi.ar
cpgram(ksi.ar$resid, main = "AR fit to residuals")
ksi.ar.resid <- as.vector(ksi.ar$resid)
plot(ksi.ar.resid, type = "l")
# ARIMA model forecasting
ksi.fore <- predict(ksi.arima, 10); ksi.fore
plot(ksi.fore$pred, type = "l")
plot(ksi.fore$se, type = "l")
ts.plot(ksi, ksi.fore$pred, ksi.fore$pred + 2*ksi.fore$se, ksi.fore$pred - 2*ksi.fore$se)
ts.plot(ksi.fore$pred, ksi.fore$pred + 2*ksi.fore$se, ksi.fore$pred - 2*ksi.fore$se)
# seasonality analysis
ksi.ts <- ts(ksi, frequency = 60); ksi.ts
ksi.stl <- stl(ksi.ts, "periodic"); ksi.stl
ksi.stl.seas <- ksi.stl$time.series[,1]  # seasonality
ksi.stl.trend <- ksi.stl$time.series[,2] # trend
ksi.stl.rem <- ksi.stl$time.series[,3]   # remainder
plot(ksi.stl.seas, type = "l")
plot(ksi.stl.trend, type = "l")
plot(ksi.stl.rem, type = "l")
# SSA analysis ofv residuals
library(Rssa)
ksi.ssa <- ssa(ksi)
ksi.ssa.series <- unlist(ksi.ssa$series)
ksi.ssa.series
plot(ksi.ssa.series)
library(kernlab)
#####################################################################
# lines 32-61 contain estimated bearing rates (sonar system output) #
#####################################################################
dBe <- SSdat[32,]
eta <- unlist(dBe)
plot(eta, type = "l")
# mean
eta.mean <- mean(eta); eta.mean
# mean root-square error
eta.mrse <- sqrt(var(eta)); eta.mrse
# analyse noise signal
# compute autocovariance
eta.acf <- acf(eta, lag.max = 3,type = "correlation")
eta.acf <- acf(eta, type = "correlation")
eta.acf <- acf(eta, type = "covariance")
eta.acf <- acf(eta, type = "partial")
eta.acf
# spectral density
eta.spec <- spectrum(eta)
# ARIMA model for residuals
eta.arima <- arima(eta, order = c(10,0,0)); eta.arima
cpgram(eta.arima$resid, main = "ARIMA fit to residuals")
# ARIMA model with Akaike's information criterion (AIC) model selection
eta.ar <- ar(eta, order.max = 10); eta.ar
cpgram(eta.ar$resid, main = "AR fit to residuals")
eta.ar.resid <- as.vector(eta.ar$resid)
plot(eta.ar.resid, type = "l")
# seasonality analysis
eta.ts <- ts(eta, frequency = 60); eta.ts
eta.stl <- stl(eta.ts, "periodic"); eta.stl
eta.stl.seas <- eta.stl$time.series[,1]  # seasonality
eta.stl.trend <- eta.stl$time.series[,2] # trend
eta.stl.rem <- eta.stl$time.series[,3]   # remainder
plot(eta.stl.seas, type = "l")
plot(eta.stl.trend, type = "l")
plot(eta.stl.rem, type = "l")
#####################################
# test for uncorrelated normal r.v. #
#####################################
psi <- rnorm(length(ksi))
plot(psi, type = "l")
# mean
psi.mean <- mean(psi); psi.mean
# mean root-square error
psi.mrse <- sqrt(var(psi)); psi.mrse
# analyse noise signal
# compute autocovariance
psi.acf <- acf(psi, type = "correlation")
psi.acf
# spectral density
psi.spec <- spectrum(psi)
# cumulative periodogram
psi.cpgram <- cpgram(psi)
#########################
# CHECK tmalib bearings #
#########################
TMAdat <- read.table("tmalib_ssda.txt", dec = ",")
TMAdat <- read.table("tmalib_ssda_1.txt", dec = ",")
TMAdat <- read.table("tmalib_ssda_2.txt", dec = ",")
Bm <- TMAdat[,1]
Be <- TMAdat[,2]
ksi <- Bm - Be
plot(ksi, type = "l")
# mean
ksi.mean <- mean(ksi); ksi.mean
# mean root-square error
ksi.mrse <- sqrt(var(ksi)); ksi.mrse
# compute autocovariance
ksi.acf <- acf(ksi, type = "correlation")
# ARIMA model for residuals
ksi.arima <- arima(ksi, order = c(10,0,0)); ksi.arima
cpgram(ksi.arima$resid, main = "ARIMA fit to residuals")
# ARIMA model with Akaike's information criterion (AIC) model selection
ksi.ar <- ar(ksi, order.max = 1); ksi.ar
cpgram(ksi.ar$resid, main = "AR fit to residuals")
# spectral density
ksi.spec <- spectrum(ksi)
# cumulative periodogram
ksi.cpgram <- cpgram(ksi)
# k-step differences (lag = k)
zeta <- diff(ksi, lag = 1, differences = 1)
plot(zeta, type = "l")
zeta.acf <- acf(zeta, type = "correlation")
zeta.cpgram <- cpgram(zeta)
zeta.ar <- ar(zeta, order.max = 30); zeta.ar
cpgram(zeta.ar$resid, main = "AR fit to residuals")
zeta.ar.resid <- as.vector(zeta.ar$resid)
plot(zeta.ar.resid, type = "l")
# seasonality analysis
ksi.ts <- ts(ksi, frequency = 60)
ksi.stl <- stl(ksi.ts, "periodic")
ksi.stl.seas <- ksi.stl$time.series[,1]  # seasonality
ksi.stl.trend <- ksi.stl$time.series[,2] # trend
ksi.stl.rem <- ksi.stl$time.series[,3]   # remainder
plot(ksi.stl.seas, type = "l")
plot(ksi.stl.trend, type = "l")
plot(ksi.stl.rem, type = "l")
#################
# SSDA AR model #
#################
SSmod <- read.table("SSDA_AR_model.txt", header = TRUE, dec = ",")
edit(SSmod)
plot(SSmod$MRSE, SSmod$C1)
plot(SSmod$BR, SSmod$C1)
plot(SSmod$MRSE, SSmod$sigma)
plot(SSmod$BR, SSmod$sigma)
plot(SSmod$BR, SSmod$MRSE)
plot(SSmod$sigma, SSmod$C1)
SSmod.C1.mean.ByBR <- tapply(SSmod$C1, SSmod$BR, mean);SSmod.C1.mean.ByBR
plot(SSmod.C1.mean.ByBR)
C1 <- mean(SSmod$C1); C1
SSmod.sigma.mean.ByBR <- tapply(SSmod$sigma, SSmod$BR, mean);SSmod.sigma.mean.ByBR
plot(SSmod.sigma.mean.ByBR)
SSmod.sigma.mean.ByMRSE <- tapply(SSmod$sigma, SSmod$MRSE, mean);SSmod.sigma.mean.ByMRSE
plot(SSmod.sigma.mean.ByMRSE)
# construct the model
SSmod.sim <- arima.sim(list(order = c(1,0,0), ar = -0.2708704), sd = 0.08, n = 200); SSmod.sim
# SSmod.sim.sd <- sqrt(var(SSmod.sim)); SSmod.sim.sd
plot(SSmod.sim, type = "l")
# check the model
# ARIMA model for residuals
SSmod.sim.arima <- arima(SSmod.sim, order = c(10,0,0)); SSmod.sim.arima
cpgram(SSmod.sim.arima$resid, main = "ARIMA fit to residuals")
# ARIMA model with Akaike's information criterion (AIC) model selection
SSmod.sim.ar <- ar(SSmod.sim, order.max = 1); SSmod.sim.ar
cpgram(SSmod.sim.ar$resid, main = "AR fit to residuals")
sd <- sqrt(0.8702); sd
