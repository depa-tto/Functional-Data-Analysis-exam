###################################################
### tecator dataset
###################################################
library("fda.usc")
data("tecator")
names(tecator)
class(tecator)

absorp <- tecator$absorp.fdata
head(absorp)
class(absorp)

head(tecator$y)
Fat20 <- ifelse(tecator$y$Fat < 20, 0, 1) * 2 + 2
absorp$names$main <- ""

dev.new( width = 150,height = 110, units = "mm")
opar <- par(mfrow = c(1 , 2))
# Figure 1 (left panel)
plot(absorp, col = Fat20)

absorp.d1 <- fdata.deriv(absorp, nderiv = 1)
# Figure 1 (right panel)
plot(absorp.d1, col = Fat20)
par(opar)

# the curves overlap, but if I compute the derivative(second graph)
# the differences can be seen more clearly

###################################################
### convert fdata class to fd class
###################################################
class(absorp.fd <- fdata2fd( absorp, type.basis = "fourier", nbasis = 15))
class(absorp.fdata <- fdata(absorp.fd)) # data discretization

dev.new( width = 150,height = 110, units = "mm")
opar <- par(mfrow = c(1 , 2))
plot(absorp.fd)
plot(absorp.fdata)
par(opar)

###################################################
### phoneme data and smoothing
###################################################
data(phoneme)
plot(phoneme$test)

learn <- phoneme$learn
l <- c(0 ,2^seq(-2, 9, len = 30))
nb <- seq(7, 31, by = 2)

out0 <- optim.basis(learn, lambda = l, numbasis = nb, type.basis = "bspline")
out0$gcv.opt
out0$numbasis.opt

opar <- par(mfrow=c(1,2))
plot(out0$fdataobj)
plot(out0$gcv["27",], type = "l")
par(opar)

out1 <- optim.np(learn, type.S = S.NW, par.CV = list(criteria = "GCV")) # non parametric karnel smoothing

out1$gcv.opt

opar <- par(mfrow=c(1,2))
plot(out1$fdataobj)
plot(out1$gcv, type = "l")
par(opar)

out2 <- optim.np(learn, type.S = S.LLR, par.CV = list(criteria = "GCV"))

out2$gcv.opt

opar <- par(mfrow=c(1,2))
plot(out2$fdataobj)
plot(out2$gcv, type = "l")
par(opar)

# Kernels
out3 <- optim.np(learn, type.S = S.KNN, h = 3:35, Ker = Ker.norm) # Normal Kernel

out4 <- optim.np(learn, type.S = S.NW, h = 3:35, Ker = Ker.tri, correl = FALSE) #Triweight Kernel

out5 <- optim.np(learn, type.S = S.NW, h = 3:35, Ker = Ker.epa, correl = FALSE) #Epanechnikov Kerne

out6 <- optim.np(learn, type.S = S.NW, h = 3:35, Ker = Ker.unif, correl = FALSE) #Uniform Kernel

###################################################
### plot GCV criteria
###################################################
dev.new(width = 150, height = 110, units = "mm")
par(mfrow = c(1,2))
contour(nb, l, out0$gcv, ylab = "Lambda", xlab = "Number of basis", 
        main = "GCV criteria by optim.basis()")
plot(out1$h, out1$gcv, type = "l", main = "GCV criteria  by optim.np() ", 
     xlab = "Bandwidth (h) values",ylab = "GCV criteria", col = 3, lwd = 2)
legend(x = 3, y = 6, legend = c("Ker.norm-S.NW", "Ker.norm-S.LLR", 
                                  "Ker.norm-S.KNN", "Ker.tri-S.NW",
                                  "Ker.epa-S.NW", "Ker.unif-S.NW"),
       box.col = "white", lwd = c(2, 2, 2), col = c(3, 4, 5, 6, 7, 8),cex = 0.75)
lines(out2$h,out2$gcv, col = 4, lwd = 2)
lines(out3$h,out3$gcv, col = 5, lwd = 2)
lines(out4$h,out4$gcv, col = 6, lwd = 2)
lines(out5$h,out5$gcv, col = 7, lwd = 2)
lines(out6$h,out6$gcv, col = 8, lwd = 2)


library(rgl)
persp3d(nb, l, out0$gcv, col="skyblue")

###################################################
### smoothing a fdata curve
###################################################
dev.new( width = 150, height = 100, units = "mm")
ind <- 11
nam <- expression( paste("Phoneme curve"[11]) )
plot(learn[ind, ], main = nam, lty = 2, lwd = 2, col = 8)
legend(x = 70, y = 19, legend = c("Curve","Bspline basis",
                                  "Ker.norm-S.NW", "Ker.norm-S.LLR", 
                                  "Ker.norm-S.KNN", "Ker.tri-S.NW",
                                  "Ker.epa-S.NW", "Ker.unif-S.NW"),
       lty = c(2, 1, 1, 1, 1, 1, 1, 1), lwd = 2, col = c(8, 1, 3, 4, 5, 6, 7, 2), box.col = "white")
lines(out0$fdata.est[ind, ], col = 1, lty = 1, lwd = 2)
lines(out1$fdata.est[ind, ], col = 3, lty = 1, lwd = 2)
lines(out2$fdata.est[ind, ], col = 4, lty = 1, lwd = 2)
lines(out3$fdata.est[ind, ], col = 5, lty = 1, lwd = 2)
lines(out4$fdata.est[ind, ], col = 6, lty = 1, lwd = 2)
lines(out5$fdata.est[ind, ], col = 7, lty = 1, lwd = 2)
lines(out6$fdata.est[ind, ], col = 2, lty = 1, lwd = 2)

# It better to choose the less noisy curves, so KNN is one of the best

###############################################################################
library(fda)
library(funData)
library(refund)

# Data
data("CanadianWeather", package = "fda")
dailyTemp <- funData(argvals = 1:365,
                     X = t(CanadianWeather$dailyAv[, , "Temperature.C"]))
monthlyPrec <- funData(argvals = 1:12,
                       X = t(CanadianWeather$monthlyPrecip))
canadWeather <- multiFunData(dailyTemp, monthlyPrec)

# Output
dailyTemp
class(dailyTemp)

# SImple plot
plot(dailyTemp, main = "Daily Temperature Data", xlab = "Day of Year",
     ylab = "Temperature in °C")

#Plot
library("ggplot2")
tempPlot <- autoplot(dailyTemp)
tempPlot + labs(title = "Daily Temperature Data",
                x = "Day of Year", y = "Temperature in °C")

# Multivariate plot
canadWeather
plot(canadWeather, obs = 26:35, lwd = 2, log = c("", "y"),
     main = c("Temperature", "Precipitation (log-scale)"),
     xlab = c("Day of Year", "Month"),
     ylab = c("Temperature in °C", "Precipitation in mm"))

# Multivariate plot
weatherPlot <- autoplot(canadWeather, obs = 26:35)

weatherPlot[[1]] <- weatherPlot[[1]] + geom_line(aes(colour = obs)) +
  labs(title = "Temperature", colour = "Weather Station",
       x = "Day of Year", y = "Temperature in °C")
weatherPlot[[2]] <- weatherPlot[[2]] + geom_line(aes(colour = obs)) +
  labs(title = "Precipitation (log-scale)", colour = "Weather Station",
       x = "Month", y = "Precipitation in mm") +
  scale_x_continuous(breaks = 1:12) +
  scale_y_log10(breaks = c(0.1, 0.5, 1, 5, 10))
gridExtra::grid.arrange(grobs = weatherPlot, nrow = 1)

# Data
install.packages('refund')
data("cd4", package = "refund")
?cd4

allArgvals <- seq(-18, 42)
argvalsList <- apply(cd4, 1, function(x) allArgvals[complete.cases(x)])
obsList <- apply(cd4, 1, function(x) x[complete.cases(x)])
cd4Counts <- irregFunData(argvals = argvalsList, X = obsList)

argvalsList
cd4Counts

# Simple plot
plot(cd4Counts, obs = 1:5, xlim = c(-18, 45), log = "y",
     main = "CD4 Counts for Individuals 1-5",
     xlab = "Month since seroconversion",
     ylab = "CD4 cell count (log-scale)")
legend("topright", legend = 1:5, col = rainbow(5), lty = 1, pch = 20,
      title = "Individual")

# Plot
cd4Plot <- autoplot(cd4Counts, obs = 1:5)
cd4Plot + geom_line(aes(colour = obs)) +
  labs(title = "CD4 Counts for Individuals 1-5", color = "Individual",
       x = "Month since seroconversion",
       y = "CD4 cell count (log-scale)") +
  scale_y_log10(breaks = seq(200, 1000, 200))

# Summary
summary(dailyTemp)
options(max.print = 24, digits = 4, scipen = 1)
summary(dailyTemp[1:6])

options(max.print = 12, digits = 7, scipen = 0)
str(cd4Counts)
summary(cd4Counts)


argvals(monthlyPrec)
names(monthlyPrec) <- names(dailyTemp)
names(monthlyPrec)
nObs(dailyTemp)
nObs(cd4Counts)
nObs(canadWeather)

nObsPoints(dailyTemp)
nObsPoints(cd4Counts)
nObsPoints(canadWeather)

dimSupp(dailyTemp)
dimSupp(cd4Counts)
dimSupp(canadWeather)

dailyTemp[1:5]
extractObs(cd4Counts, obs = 1:8, argvals = -18:0)

# Coercion 
as.irregFunData(dailyTemp)
as.data.frame(cd4Counts)

# Mathematical operations
op1 <- 9/5 * dailyTemp + 32
opar <- par(mfrow = c(1,2))
plot(dailyTemp)
plot(op1)
par(opar)

op2 <- log(cd4Counts)
opar <- par(mfrow = c(1,2))
plot(cd4Counts)
plot(op2)
par(opar)

op3 <- canadWeather - meanFunction(canadWeather)
plot(canadWeather)
plot(op3)


tensorData <- tensorProduct(dailyTemp, monthlyPrec)
tensorData
plot(tensorData, obs = 1)
plot(tensorData, obs = 2)

# Simualtion study
simUniv1D <- simFunData(N = 8, argvals = seq(0, 1, 0.01),
                        eFunType = "Fourier", eValType = "linear", M = 10)
simUniv1D
class(simUniv1D)

plot(simUniv1D$simData)

argvalsList <- list(seq(0, 1, 0.01), seq(-0.5, 0.5, 0.01))
simUniv2D <- simFunData(N = 5, argvals = argvalsList,
                        eFunType = c("Wiener", "Fourier"), eValType = "linear", M = c(10, 12))

plot(simUniv2D$simData, obs = 1)
plot(simUniv2D$simData, obs = 2)
plot(simUniv2D$simData, obs = 3)
plot(simUniv2D$simData, obs = 4)
plot(simUniv2D$simData, obs = 5)

#Calling simMultiFunData with the option "split" constructs multivariate eigenfunctions by
#splitting orthonormal functions into p pieces and shifting them to where the elements should
#be defined. 
argvalsList <- list(seq(-0.5, 0.5, 0.01), seq(0, 1, 0.01))
simMultiSplit <- simMultiFunData(N = 7, argvals = argvalsList,
                                 eFunType = "Fourier", eValType = "linear", M = 10, type = "split")

# As an alternative, multivariate eigenfunctions can be constructed as weighted versions of
# univariate eigenfunctions. 
argvalsList <- list(list(seq(-0.5, 0.5, 0.01)), list(seq(0, 1, 0.01), seq(-1, 1, 0.01)))
simMultiWeight <- simMultiFunData(N = 5, argvals = argvalsList,
                                  eFunType = list("Fourier", c("Wiener", "Poly")),
                                  eValType = "exponential", M = list(12, c(4, 3)), type = "weighted")

sim1Derr <- addError(simUniv1D$simData, sd = 0.5)
sim1Dsp <- sparsify(simUniv1D$simData, minObs = 5, maxObs = 10)

plot(sim1Derr)
plot(simUniv1D$simData)
points(sim1Dsp)

sim2Derr <- addError(simMultiWeight$simData, sd = c(0.5, 0.3))
sparsify(simMultiSplit$simData, minObs = c(5, 50), maxObs = c(10, 80))



plot(sim2Derr, obs = 1)
plot(sim2Derr, obs = 2)
plot(sim2Derr, obs = 3)
plot(sim2Derr, obs = 4)
plot(sim2Derr, obs = 5)













