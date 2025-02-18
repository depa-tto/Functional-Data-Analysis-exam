### Source: https://rdrr.io/rforge/fda/src/inst/scripts/fdarm-ch04.R

###
### Ramsay, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###

#  Remarks and disclaimers

#  These R commands are either those in this book, or designed to
#  otherwise illustrate how R can be used in the analysis of functional
#  data.
#  We do not claim to reproduce the results in the book exactly by these
#  commands for various reasons, including:
#    -- the analyses used to produce the book may not have been
#       entirely correct, possibly due to coding and accuracy issues
#       in the functions themselves
#    -- we may have changed our minds about how these analyses should be
#       done since, and we want to suggest better ways
#    -- the R language changes with each release of the base system, and
#       certainly the functional data analysis functions change as well
#    -- we might choose to offer new analyses from time to time by
#       augmenting those in the book
#    -- many illustrations in the book were produced using Matlab, which
#       inevitably can imply slightly different results and graphical
#       displays
#    -- we may have changed our minds about variable names.  For example,
#       we now prefer "yearRng" to "yearRng" for the weather data.
#    -- three of us wrote the book, and the person preparing these scripts
#       might not be the person who wrote the text
#  Moreover, we expect to augment and modify these command scripts from time
#  to time as we get new data illustrating new things, add functionality
#  to the package, or just for fun.

# Some parts of the script are changed by Jurgita Markevi i t  (2021)

###
### ch. 4  How to Build Functional Data Objects
###

#  load the fda package

library(fda)

#  display the data files associated with the fda package

data(package='fda')

#  start the HTML help system if you are connected to the Internet, in
#  order to open the R-Project documentation index page in order to obtain
#  information about R or the fda package.

help.start()

##
## Section 4.1 Adding Coefficients to Bases to Define Functions
##

#  4.1.1 Coefficient Vectors, Matrices and Arrays

daybasis65 = create.fourier.basis(c(0,365), 65)
# dummy coefmat
coefmat = matrix(0, 65, 35, dimnames=list(
  daybasis65$names, CanadianWeather$place) )
tempfd. = fd(coefmat, daybasis65)
plot(tempfd.)

# 4.1.2 Labels for Functional Data Objects

fdnames      = list("Age (years)", "Child", "Height (cm)")

# or

fdnames      = vector('list', 3)
fdnames[[1]] = "Age (years)"
fdnames[[2]] = "Child"
fdnames[[3]] = "Height (cm)"

station      = vector('list', 35)
station[[1]]= "St. Johns"
#.
#.
#.
station[[35]] = "Resolute"

# Or:

station = as.list(CanadianWeather$place)

fdnames = list("Day", "Weather Station" = station,
               "Mean temperature (deg C)")

##
## 4.2 Methods for Functional Data Objects
##

#  Two order 2 splines over unit interval

unitRng = c(0,1)
bspl2 = create.bspline.basis(unitRng, norder=2)
plot(bspl2, lwd=2)

#  a pair of straight lines

tstFn1 = fd(c(-1, 2), bspl2)
tstFn2 = fd(c(1, 3), bspl2)


opar <- par(mfrow=c(3,2))

plot(tstFn1,   lwd=2, xlab="", ylab="Line 1")
plot(tstFn2,   lwd=2, xlab="", ylab="Line 2")

#  sum of these straight lines
fdsumobj = tstFn1+tstFn2
plot(fdsumobj, lwd=2, xlab="", ylab="Line1 + Line 2")

#  difference between these lines
fddifobj = tstFn2-tstFn1
plot(fddifobj, lwd=2, xlab="", ylab="Line2 - Line 1")

# multiplication of two lines
fdprdobj = tstFn1 * tstFn2
plot(fdprdobj, lwd=2, xlab="", ylab="Line2 * Line 1")

#  square of a straight line
fdsqrobj = tstFn1^2
plot(fdsqrobj, lwd=2, xlab="", ylab="Line 1 ^2")
par(opar)

#devide two lines
fdprdobj = tstFn1 / tstFn2
# it is not allowed

tt <- seq(0,1,by=0.05)
ev1 <- eval.fd(tt, tstFn1)
ev2 <- eval.fd(tt, tstFn2)
div <- ev1/ev2

#  square root of a line with negative values:  illegal

a = 0.5
fdrootobj = tstFn1^a
#Error in `^.fd`(tstFn1, a) :
#  There are negative values and the power is a positive fraction.

#  square root of a square:  this illustrates the hazards of
#  fractional powers when values are near zero.  The right answer is
#  two straight line segments with a discontinuity in the first
#  derivative.  It would be better to use order two splines and
#  put a knot at the point of discontinuity, but the power method
#  doesn't know how to do this.

opar <- par(mfrow=c(3,1))
fdrootobj = fdsqrobj^a
plot(tstFn1,    lwd=2, xlab="", ylab="Line 1")
plot(tstFn2,    lwd=2, xlab="", ylab="Line 2")
plot(fdrootobj, lwd=2, xlab="", ylab="sqrt(fdsqrobj)")
par(opar)

#  square root of a quadratic without values near zero:  no problem

fdrootobj = (fdsqrobj + 1)^a
opar <- par(mfrow=c(2,1))
plot(fdsqrobj + 1,    lwd=2, xlab="", ylab="fdsqrobj + 1")
plot(fdrootobj, lwd=2, xlab="", ylab="sqrt(fdsqrobj + 1)")
par(opar)

#  reciprocal of a function with zero values:  illegal operation

a    = (-1)
fdinvobj = tstFn1^a
# Error in `^.fd`(tstFn1, a) :
#   There are zero or negative values and the power is negative.

#  reciprocal of a function with near zero values:  a foolish thing
#  to do and the power function fails miserably

opar <- par(mfrow=c(2,1))
fdinvobj = fdsqrobj^a
plot(fdsqrobj, lwd=2, xlab="", ylab="fdsqrobj")
plot(fdinvobj, lwd=2, xlab="", ylab="1/fdsqrobj")
par(opar)

#  reciprocal of a positive function with no values near zero

opar <- par(mfrow=c(2,1))
fdinvobj = (fdsqrobj+1)^a
plot(fdsqrobj + 1, lwd=2, xlab="", ylab="fdsqrobj + 1")
plot(fdinvobj,     lwd=2, xlab="", ylab="1/(fdsqrobj+1)")
par(opar)

#  near reciprocal of a positive function with no values near zero

a = -0.99
opar <- par(mfrow=c(2,1))
fdpowobj = (fdsqrobj+1)^a
plot(fdsqrobj + 1, lwd=2, xlab="", ylab="fdsqrobj + 1")
plot(fdpowobj,     lwd=2, xlab="", ylab="(fdsqrobj+1)^(-0.99)")
par(opar)

#
# compute mean temperature in two ways and plot the difference
#

yearRng = c(0,365)

Tempbasis = create.fourier.basis(yearRng, 65)
Tempfd = smooth.basis(day.5,
                      CanadianWeather$dailyAv[,,'Temperature.C'], Tempbasis)$fd
meanTempfd = mean.fd(Tempfd)

plot(Tempfd)
lines(meanTempfd, col = "red", lwd = 2)

sumTempfd  = sum(Tempfd)

plot((meanTempfd-sumTempfd*(1/35)))

# round off error, as it should be.

#  plot the temperature for Resolute and add the Canadian mean

plot(Tempfd[35], lwd=1, ylim=c(-35,20))
lines(meanTempfd, lty=2)

#  evaluate the derivative of mean temperature and plot

DmeanTempVec = eval.fd(day.5, meanTempfd, 1)
plot(day.5, DmeanTempVec, type='l')

#  evaluate and plot the harmonic acceleration of mean temperature

harmaccelLfd = vec2Lfd(c(0,c(2*pi/365)^2, 0), c(0, 365))
LmeanTempVec = eval.fd(day.5, meanTempfd, harmaccelLfd)

par(mfrow=c(1,1))
plot(day.5, LmeanTempVec, type="l", cex=1.2,
     xlab="Day", ylab="Harmonic Acceleration")
abline(h=0)

#  plot Figure 4.1

dayOfYearShifted = c(182:365, 1:181)

tempmat   = daily$tempav[dayOfYearShifted, ]
tempbasis = create.fourier.basis(yearRng,65)

temp.fd = smooth.basis(day.5, tempmat, tempbasis)$fd

temp.fd$fdnames = list("Day (July 2 to June 30)",
                       "Weather Station",
                       "Mean temperature (deg. C)")

plot(temp.fd, lwd=2, xlab='Day (July 1 to June 30)',
     ylab='Mean temperature (deg. C)')

#
# Section 4.2.1 Illustration: Sinusoidal Coefficients
#

# Figure 4.2

basis13  = create.bspline.basis(c(0,10), 13)
tvec     = seq(0,1,len=13)
sinecoef = sin(2*pi*tvec)
sinefd   = fd(sinecoef, basis13, list("t","","f(t)"))
op       = par(cex=1.2)
plot(sinefd, lwd=2)
points(tvec*10, sinecoef, lwd=2)
par(op)

##
## Section 4.3 Smoothing using Regression Analysis
##

# Section 4.3.1 Plotting the January Thaw

# Figure 4.3

# This assumes the data are in "MtlDaily.txt"
# in the working directory getwd();
# first create it and put it there
cat(MontrealTemp, file='MtlDaily.txt')

MtlDaily = matrix(scan("MtlDaily.txt",0),34,365)
thawdata = t(MtlDaily[,16:47])

daytime  = ((16:47)+0.5)
plot(daytime, apply(thawdata,1,mean), "b", lwd=2,
     xlab="Day", ylab="Temperature (deg C)", cex=1.2)

# Figure 4.4
# Now we can compute coefficients for our functional data object by the usual equations
# for regression coefficients, b = (X'X)^(-1) X'y; and construct a functional data
# object by combining them with our basis object. 
# Wwe do see a fair number of them peaking between January 20 and 25, and a 
# few others with later peaks as well.

thawbasis    = create.bspline.basis(c(16,48),7)
thawbasismat = eval.basis(thawbasis, daytime)

thawcoef = solve(crossprod(thawbasismat),
                 crossprod(thawbasismat,thawdata))
thawfd   = fd(thawcoef, thawbasis,
              list("Day", "Year", "Temperature (deg C)"))
plot(thawfd, lty=1, lwd=2, col=1)

# Figure 4.5

plotfit.fd(thawdata[,1], daytime, thawfd[1],
           lty=1, lwd=2, main='')

plotfit.fd(thawdata, daytime, thawfd,
           lty=1, lwd=2, main='')

##
## Section 4.4 The Linear Differential Operator or Lfd Class
##

#---
  #http://faculty.bscb.cornell.edu/~hooker/FDA2008/Lecture7_refinery.R

# The first thing we need to do is to define a linear differential operator 
# (Lfd) object. This is a list of functional data objects defining the terms 
# on the right hand side of
#
#  D^m x + b_{m-1}(t) D^{m-1} x + ... + b_0(t) x

?Lfd

# The simplest Lfds just penalize the derivative x, and we can define
# them by 

D2lfd = int2Lfd(2)
D2lfd

# This says that b_0 = b_1 = 0, which we can see by looking internally at
# D2lfd we see

names(D2lfd)

plot(D2lfd$bwtlist[[2]])

# Instead of evaluating a derivative, we can always evaluate an Lfd
# applied to a function

plot(temp.fd,Lfdobj=D2lfd)
plot(temp.fd,Lfdobj=1)

# Now we need to do some smoothing. 
# In order to set up the smoothing operator, we need to define an 
# fdPar (functional parameter) object. This is just a list of various
# quantities that avoids the need for many argument values

ageRng  = c(1,18)
age     = growth$age
agefine = seq(1,18,len=501)

#  set up order 6 spline basis with 12 basis functions for
#  fitting the growth data so as to estimate acceleration

nbasis = 12;
norder =  6;
heightbasis12 = create.bspline.basis(ageRng, nbasis, norder)

#  fit the data by least squares

basismat   = eval.basis(age, heightbasis12)
heightmat  = growth$hgtf
heightcoef = lsfit(basismat, heightmat, intercept=FALSE)$coef

D2fdPar = fdPar(heightbasis12,Lfdobj=int2Lfd(2),lambda=1e4)

# The functional parameter object holds a basis, a Lfd and a value for lambda. 
# We can now make a new call to 

syfd = smooth.basis(age, heightmat, D2fdPar)

plotfit.fd(heightmat,age,syfd$fd)


# Now we can examine how this changes with lambda. It is most useful to
# vary lambda on the logarithmic scale. Note that large values of 
# lambda mean more smoothing

# We'll also keep track of gcv, df and sse

gcv = rep(0,21)
df = rep(0,21)
sse = rep(0,21)


for(i in 1:21){
  lambda=10^{i-10}
  tD2fdPar = fdPar(heightbasis12,Lfdobj=int2Lfd(2),lambda=lambda)
  
  tyfd = smooth.basis(age, heightmat, tD2fdPar)
  
  gcv[i] = sum(tyfd$gcv)
  df[i] = tyfd$df
  sse[i] = tyfd$SSE
}

# And we'll plot some results

plot(-10:0,df[1:11],type='l',xlab='log lambda',ylab='df',cex.lab=1.5)
plot(-10:0,sse[1:11],type='l',xlab='log lambda',ylab='sse',cex.lab=1.5)
plot(-10:0,gcv[1:11],type='l',xlab='log lambda',ylab='gcv',cex.lab=1.5)


# We'll try the same thing with very small lambda

D2fdPar1 = fdPar(heightbasis12,Lfdobj=int2Lfd(2),lambda=1e-6)
syfd1 = smooth.basis(age,heightmat,D2fdPar1)
plotfit.fd(heightmat,age,syfd1$fd)

# Or very large lambda

D2fdPar2 = fdPar(heightbasis12,Lfdobj=int2Lfd(2),lambda=1e10)
syfd2 = smooth.basis(age,heightmat,D2fdPar2)
plotfit.fd(heightmat,age,syfd2$fd)

# Or min GCV

D2fdPar3 = fdPar(heightbasis12,Lfdobj=int2Lfd(2),lambda=0.001)
syfd3 = smooth.basis(age,heightmat,D2fdPar3)
plotfit.fd(heightmat,age,syfd3$fd)
#---


omega           = 2*pi/365
thawconst.basis = create.constant.basis(thawbasis$rangeval)

betalist       = vector("list", 3)
betalist[[1]]  = fd(0, thawconst.basis)
betalist[[2]]  = fd(omega^2, thawconst.basis)
betalist[[3]]  = fd(0, thawconst.basis)
harmaccelLfd.  = Lfd(3, betalist)

accelLfd = int2Lfd(2)

harmaccelLfd.thaw = vec2Lfd(c(0,omega^2,0), thawbasis$rangeval)
all.equal(harmaccelLfd.[-1], harmaccelLfd.thaw[-1])

class(accelLfd)
class(harmaccelLfd)

Ltempmat  = eval.fd(day.5, temp.fd, harmaccelLfd)

D2tempfd = deriv.fd(temp.fd, 2)
Ltempfd  = deriv.fd(temp.fd, harmaccelLfd)

##
## Section 4.5 Bivariate Functional Data Objects:
##             Functions of Two Arguments
##

Bspl2 = create.bspline.basis(nbasis=2, norder=1)
Bspl3 = create.bspline.basis(nbasis=3, norder=2)

corrmat  = array(1:6/6, dim=2:3)
bBspl2.3 = bifd(corrmat, Bspl2, Bspl3)

### Source : https://rdrr.io/rforge/fda/src/inst/scripts/fdarm-ch05.R

###
### Ramsay, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###

#  Remarks and disclaimers

#  These R commands are either those in this book, or designed to
#  otherwise illustrate how R can be used in the analysis of functional
#  data.
#  We do not claim to reproduce the results in the book exactly by these
#  commands for various reasons, including:
#    -- the analyses used to produce the book may not have been
#       entirely correct, possibly due to coding and accuracy issues
#       in the functions themselves
#    -- we may have changed our minds about how these analyses should be
#       done since, and we want to suggest better ways
#    -- the R language changes with each release of the base system, and
#       certainly the functional data analysis functions change as well
#    -- we might choose to offer new analyses from time to time by
#       augmenting those in the book
#    -- many illustrations in the book were produced using Matlab, which
#       inevitably can imply slightly different results and graphical
#       displays
#    -- we may have changed our minds about variable names.  For example,
#       we now prefer "yearRng" to "yearRng" for the weather data.
#    -- three of us wrote the book, and the person preparing these scripts
#       might not be the person who wrote the text
#  Moreover, we expect to augment and modify these command scripts from time
#  to time as we get new data illustrating new things, add functionality
#  to the package, or just for fun.

# Some parts of this script corrected by Jurgita Markeviciute (2021)

###
### ch. 5.  Smoothing: Computing Curves from Noisy Data
###

#  load the fda package

library(fda)

###############################################################################
##
## Section 5.1.  Regression Splines: Smoothing by Regression Analysis
##

#  -------------------  Smoothing the growth data  ------------------------

#  define the range of the ages and set up a fine mesh of ages

rm(list=ls())

ageRng  = c(1,18)
age     = growth$age
agefine = seq(1,18,len=501)

#  set up order 6 spline basis with 12 basis functions for
#  fitting the growth data so as to estimate acceleration

nbasis = 12;
norder =  6;
heightbasis12 = create.bspline.basis(ageRng, nbasis, norder)

#  fit the data by least squares

basismat   = eval.basis(age, heightbasis12)
heightmat  = growth$hgtf
heightcoef = lsfit(basismat, heightmat, intercept=FALSE)$coef

fdnames      = vector('list', 3)
fdnames[[1]] = "Age (years)"
fdnames[[2]] = "Child"
fdnames[[3]] = "Height (cm)"

hgtfd <- fd(heightcoef, heightbasis12, fdnames)
plot(hgtfd)
str(hgtfd)
names(hgtfd)

#  fit the data using function smooth_basis, which does the same thing.

heightList = smooth.basis(age, heightmat, heightbasis12)
heightfd   = heightList$fd
plot(heightfd)
str(heightfd)
names(heightfd)

height.df  = heightList$df
height.df

height.gcv = heightList$gcv
height.gcv

heightbasismat = eval.basis(age, heightbasis12)
y2cMap         = solve(crossprod(heightbasismat)) %*% t(heightbasismat)

coefest <- y2cMap %*% heightmat
dcoef <- coefest - heightcoef
dcoef

# Hat matrix H
# Phi (Phi^T Phi)^{-1} Phi^T
# with penalty: Phi (Phi^T Phi + lambda R)^{-1} Phi^T
hatH = heightbasismat %*% 
  solve(crossprod(heightbasismat)) %*% t(heightbasismat)


##
## Section 5.2.  Data Smoothing with Roughness Penalties
##

# section 5.2.2 The Roughness Penalty Matrix R

#  ---------------  Smoothing the Canadian weather data  ------------------

#  define a Fourier basis for daily temperature data

yearRng   = c(0,365)
nbasis    = 65
tempbasis = create.fourier.basis(yearRng,nbasis)

#  define the harmonic acceleration operator

#  When the coefficients of the linear differential operator
#  are constant, the Lfd object can be set up more simply
#  by using function vec2Lfd as follows:

#  The first argument is a vector of coefficients for the
#  operator, and the second argument is the range over which
#  the operator is defined.

harmaccelLfd = vec2Lfd(c(0,(2*pi/365)^2,0), yearRng)

#Lx(t) = b_0(t) x(t) + b_1(t)Dx(t) + ... + b_{m-1}(t) D^{m-1} x(t) + D^mx(t)
constb <- create.constant.basis(yearRng)
bwlist <- vector("list", 3)
for (j in 1:3) bwlist[[j]] <- fd(c(0,(2*pi/365)^2,0)[j], constb)
bwlist
lfdobj <- Lfd(3, bwlist)

#  compute the penalty matrix R
# R = integral phi(t) phi^T(t) dt
# Estimator = (Phi^T Phi + lanbda R)^{-1} Phi^T y

Rmat = eval.penalty(tempbasis, harmaccelLfd)

# section 5.2.4 Defining Smoothing by Functional Parameter Objects

#  -------  Smoothing the growth data with a roughness penalty  -----------

#  set up a basis for the growth data
#  with knots at ages of height measurement

norder      = 6
nbasis      = length(age) + norder - 2
heightbasis = create.bspline.basis(ageRng, nbasis, norder, age)

#  define a functional parameter object for smoothing

heightLfd    = 4
heightlambda = 0.01
heightfdPar  = fdPar(heightbasis, heightLfd, heightlambda)

#  smooth the data

heightfdSmooth = smooth.basis(age, heightmat, heightfdPar)
heightfdpen       = heightfdSmooth$fd

opar <- par(mfrow=c(3,2))
plot(heightfd)
plot(heightfdpen)
plot(heightfd, Lfdobj = 1)
plot(heightfdpen, Lfdobj = 1)
plot(heightfd, Lfdobj = 2)
plot(heightfdpen, Lfdobj = 2)
par(opar)


# section 5.2.5 Choosing Smoothing Parameter lambda

loglam         = seq(-6, 0, 0.25)
Gcvsave        = rep(NA, length(loglam))
names(Gcvsave) = loglam
Dfsave         = Gcvsave
for(i in 1:length(loglam)){
  hgtfdPari  = fdPar(heightbasis, Lfdobj=4, 10^loglam[i])
  hgtSm.i    = smooth.basis(age, heightmat, hgtfdPari)
  Gcvsave[i] = sum(hgtSm.i$gcv)
  Dfsave[i]  = hgtSm.i$df
}

# Figure 5.1.

par(mfrow=c(1,1))
plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
     ylab=expression(GCV(lambda)), lwd=2 )
abline(h=min(Gcvsave), col=2, lty=2)
abline(v=loglam[which.min(Gcvsave)], col=2, lty=2)
abline(v=loglam[which.min(Gcvsave)+1], col=4, lty=2)
abline(v=loglam[which.min(Gcvsave)-2], col=4, lty=2)

##
## 5.3.  Case Study: The Log Precipitation Data
##

#  organize data to have winter in the center of the plot

dayOfYearShifted = c(182:365, 1:181)

dim(CanadianWeather$dailyAv)

logprecav = CanadianWeather$dailyAv[dayOfYearShifted, , 'log10precip']

#  set up a saturated basis: as many basis functions as observations

nbasis   = 365
daybasis = create.fourier.basis(yearRng, nbasis)

#  set up the harmonic acceleration operator

Lcoef        = c(0,(2*pi/diff(yearRng))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, yearRng)

#  step through values of log(lambda)

loglam        = seq(4,9,0.25)
nlam          = length(loglam)
dfsave        = rep(NA,nlam)
names(dfsave) = loglam
gcvsave       = dfsave
for (ilam in 1:nlam) {
  cat(paste('log10 lambda =',loglam[ilam],'\n'))
  lambda        = 10^loglam[ilam]
  fdParobj      = fdPar(daybasis, harmaccelLfd, lambda)
  smoothlist    = smooth.basis(day.5, logprecav,
                               fdParobj)
  dfsave[ilam]  = smoothlist$df
  gcvsave[ilam] = sum(smoothlist$gcv)
}

# Figure 5.2.

plot(loglam, gcvsave, type='b', lwd=2, ylab='GCV Criterion',
     xlab=expression(log[10](lambda)) )

#  smooth data with minimizing value of lambda

lambda      = 1e6
fdParobj    = fdPar(daybasis, harmaccelLfd, lambda)
logprec.fit = smooth.basis(day.5, logprecav, fdParobj)
logprec.fd  = logprec.fit$fd
fdnames     = list("Day (July 1 to June 30)",
                   "Weather Station" = CanadianWeather$place,
                   "Log 10 Precipitation (mm)")
logprec.fd$fdnames = fdnames

#  plot the functional data object

plot(logprec.fd, lwd=2)

fulfd <- smooth.basis(day.5, logprecav, daybasis)

opar <- par(mfrow=c(2,1))
plot(logprec.fd, lwd=2)
plot(fulfd$fd, lwd=2)
par(opar)


# Compare with actual values

plotfit.fd(logprecav, day.5, logprec.fd)

# plotfit.fd:  Pauses between plots
# *** --->>> input required (e.g., click on the plot)
#            to advance to the next plot

##
## Section 5.4 Positive, Monotone, Density
##             and Other Constrained Functions
##

#   ----------------  Positive smoothing of precipitation  ----------------

lambda      = 1e3
WfdParobj   = fdPar(daybasis, harmaccelLfd, lambda)
VanPrec     = CanadianWeather$dailyAv[
  dayOfYearShifted, 'Vancouver', 'Precipitation.mm']
VanPrecPos  = smooth.pos(day.5, VanPrec, WfdParobj)
Wfd         = VanPrecPos$Wfdobj
Wfd$fdnames = list("Day (July 1 to June 30)",
                   "Weather Station" = CanadianWeather$place,
                   "Log 10 Precipitation (mm)")

precfit = exp(eval.fd(day.5, Wfd))

plot(day.5, VanPrec, type="p", cex=1.2,
     xlab="Day (July 1 to June 30)",
     ylab="Millimeters",
     main="Vancouver's Precipitation")
lines(day.5, precfit,lwd=2)

#  ------------  5.4.2.1  Monotone smoothing  of the tibia data  -----------

#  set up the data for analysis

day    = infantGrowth[, 'day']
tib    = infantGrowth[, 'tibiaLength']
n      = length(tib)

#  a basis for monotone smoothing

nbasis = 42
Wbasis   = create.bspline.basis(c(1,n), nbasis)

#  the fdPar object for smoothing

Wfd0     = fd(matrix(0,nbasis,1), Wbasis)
WfdPar   = fdPar(Wfd0, 2, 1e-4)

#  smooth the data

result   = smooth.monotone(day, tib, WfdPar)
Wfd      = result$Wfd
beta     = result$beta

#  compute fit and derivatives of fit

dayfine  = seq(1,n,len=151)
tibhat   = beta[1]+beta[2]*eval.monfd(dayfine ,Wfd)
Dtibhat  =        beta[2]*eval.monfd(dayfine, Wfd, 1)
D2tibhat =        beta[2]*eval.monfd(dayfine, Wfd, 2)

#  plot height

op = par(mfrow=c(3,1), mar=c(5,5,3,2), lwd=2)

plot(day, tib, type = "p", cex=1.2, las=1,
     xlab="Day", ylab='', main="Tibia Length (mm)")
lines(dayfine, tibhat, lwd=2)

#  plot velocity

plot(dayfine, Dtibhat, type = "l", cex=1.2, las=1,
     xlab="Day", ylab='', main="Tibia Velocity (mm/day)")

#  plot acceleration

plot(dayfine, D2tibhat, type = "l", cex=1.2, las=1,
     xlab="Day", ylab='', main="Tibia Acceleration (mm/day/day)")
lines(c(1,n),c(0,0),lty=2)

par(op)

# ---------  5.4.2.2  Monotone smoothing the Berkeley female data  --------

##
##  Compute the monotone smoothing of the Berkeley female growth data.
##

#  set up ages of measurement and an age mesh

age     = growth$age
nage    = length(age)
ageRng  = range(age)
nfine   = 101
agefine = seq(ageRng[1], ageRng[2], length=nfine)

#  the data

hgtf   = growth$hgtf
ncasef = dim(hgtf)[2]

#  an order 6 bspline basis with knots at ages of measurement

norder = 6
nbasis = nage + norder - 2
wbasis = create.bspline.basis(ageRng, nbasis, norder, age)

#  define the roughness penalty for function W

Lfdobj    = 3          #  penalize curvature of acceleration
lambda    = 10^(-0.5)  #  smoothing parameter
cvecf     = matrix(0, nbasis, ncasef)
Wfd0      = fd(cvecf, wbasis)
growfdPar = fdPar(Wfd0, Lfdobj, lambda)

#  monotone smoothing

growthMon = smooth.monotone(age, hgtf, growfdPar)

# (wait for an iterative fit to each of 54 girls)

Wfd        = growthMon$Wfd
betaf      = growthMon$beta
hgtfhatfd  = growthMon$yhatfd

#  Set up functional data objects for the acceleration curves
#  and their mean.  Suffix UN means "unregistered".

accelfdUN     = deriv.fd(hgtfhatfd, 2)
accelmeanfdUN = mean.fd(accelfdUN)

#  plot unregistered curves

par(ask=FALSE)
plot(accelfdUN, xlim=ageRng, ylim=c(-4,3), lty=1, lwd=2,
     cex=2, xlab="Age", ylab="Acceleration (cm/yr/yr)")

