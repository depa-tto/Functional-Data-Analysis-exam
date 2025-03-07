### Source: https://rdrr.io/rforge/fda/src/inst/scripts/fdarm-ch03.R

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
### ch. 3.  How to specify basis systems for building functions
###

#  load the fda package

library(fda)

#  display the data files associated with the fda package

data(package='fda')

##
## Section 3.1 Basis Function Systems for Constructing Functions
##

unitRng = c(0,1)
unitRng

const.basis   = create.constant.basis(unitRng)
plot(const.basis)

monom.basis   = create.monomial.basis(unitRng, nbasis=5)
plot(monom.basis)

fourier.basis = create.fourier.basis(unitRng, nbasis=5, period=1)
plot(fourier.basis)
fourier.basis = create.fourier.basis(unitRng, nbasis=5, period=0.5)
plot(fourier.basis)
fourier.basis = create.fourier.basis(unitRng, nbasis=5, period=2)
plot(fourier.basis)

bspline.basis = create.bspline.basis(unitRng, nbasis=5, 
                                     norder=2, 
                                     breaks=seq(0, 1, length=5) )
plot(bspline.basis)

##
## Section 3.2 Fourier Series for Periodic Data and Functions
##

yearRng = c(0,365)

daybasis65 = create.fourier.basis(yearRng, 65)
plot(daybasis65)

tt <- seq(0, 365, by = 2)
fvals <- eval.basis(tt, daybasis65)
matplot(tt,fvals,type='l')

plot(tt,fvals[,40],type='l',lwd=1.5, lty=1)
lines(tt,fvals[,50],col=3,lwd=1.5, lty=2)
lines(tt,fvals[,30],col=4,lwd=1.5, lty=3)
lines(tt,fvals[,60],col=5,lwd=1.5, lty=4)
lines(tt,fvals[,10],col=6,lwd=1.5, lty=5)

T_size        = 500
daybasis.T = create.fourier.basis(yearRng, 3, period=T_size)
plot(daybasis.T)

# dropindthat contains a vector of indices 
# of basis functions to remove from the final series
zerobasis  = create.fourier.basis(yearRng, 65, dropind=1)
all.equal(zerobasis, daybasis65)
str(zerobasis)
str(daybasis65)

sinbasis <- create.fourier.basis(yearRng, 65, dropind=seq(3,65,by=2))
cosbasis <- create.fourier.basis(yearRng, 65, dropind=seq(2,65,by=2))
str(sinbasis)
str(cosbasis)
plot(sinbasis)
plot(cosbasis)

opar <- par(mfrow=c(2,1))
plot(zerobasis)
plot(daybasis65)
par(opar)


f.ex1 <- create.fourier.basis(unitRng, 5)
f.ex2 <- create.fourier.basis(unitRng, nbasis=5, dropind = 1)

opar <- par(mfrow=c(2,1))
plot(f.ex1)
plot(f.ex2, col=2:5, lty=2:5)
par(opar)


help(create.fourier.basis)

##
## Section 3.3 Spline series for Non-periodic Data and Functions
##

#  section 3.3.3 Examples

#  order 4 spline, one interior knot

bspline4 = create.bspline.basis(breaks=c(0, .5, 1))

# knots and knots.fd do not work for the moment
knots.fd(bspline4, interior=FALSE)
?knots.fd

# number of basis functions = order + number of interior knots.
# if we define a function over [0,1] with a single interior break point at 0.5 with cubic 
# splinebasis (order 4), then the knots are (0, 0, 0, 0, 0.5, 1, 1, 1, 1)
nord <- norder(bspline4)
rng <- bspline4$rangeval 
int <- bspline4$params
nord; rng; int
allKnots <- c(rep(rng[1], nord), int, rep(rng[2], nord))
allKnots

# we may write our own function for knots:
knots.fda <- function(Fn, interior=TRUE) {
  if(!class(bspline4) == "basisfd") stop("Object mus be of class basisfd")
  
  int <- Fn$params
  if(interior == TRUE) return(int)
  
  if(interior != TRUE) {
    nord <- norder(Fn)
    rng <- Fn$rangeval
    allKnots <- c(rep(rng[1], nord), int, rep(rng[2], nord))
    return(allKnots)
  }
}

knots.fda(bspline4, interior = TRUE)
knots.fda(bspline4, interior = FALSE)

plot(bspline4, lwd=2)

#  order 2 spline, one interiot knot

bspline2 = create.bspline.basis(breaks=c(0, .5, 1), norder=2)
knots.fda(bspline2, interior=FALSE)

plot(bspline2, lwd=2)

#  order 2 spline,  2 equal interior knots

bspline2.2 = create.bspline.basis(breaks=c(0, .5, .5, 1), norder=2)
knots.fda(bspline2.2, interior=FALSE)

plot(bspline2.2, lwd=2)

#  order 4 spline, 3 equal interior knots

bspline4.2 = create.bspline.basis(breaks=c(0, .5, .5, .5, 1), norder=2)
knots.fda(bspline4.3, interior=FALSE)

plot(bspline4.3, lwd=2)

#  section 3.3.4 B-Splines

splinebasis = create.bspline.basis(c(0,10), 13)
norder(splinebasis)
knots.fda(splinebasis)

# Figure 3.1

plot(splinebasis, xlab='t', ylab='Bspline basis functions B(t)', 
     las=1, lwd=2)

# Figure 3.2

#creat basis functions
basis2 = create.bspline.basis(c(0,2*pi), 5, 2)
basis3 = create.bspline.basis(c(0,2*pi), 6, 3)
basis4 = create.bspline.basis(c(0,2*pi), 7, 4)
plot(basis2)
plot(basis4)

#simulate data
theta     = seq(0, 2*pi, length=201)
sin.theta = sin(theta)
plot(theta, sin.theta, type="l")

#create functional data from simulated data and basis functions
sin2 = Data2fd(theta, sin.theta, basis2)
sin3 = Data2fd(theta, sin.theta, basis3)
sin4 = Data2fd(theta, sin.theta, basis4)
plot(sin2)

#evalute function values at certain points
sin2.theta = predict(sin2, theta)
sin2.theta2 <- eval.fd(theta, sin2)
head(cbind(sin2.theta, sin2.theta2))

sin3.theta = predict(sin3, theta)
sin4.theta = predict(sin4, theta)


sinRng = range(sin2.theta)
pi3    = ((1:3)*pi/2)

op = par(mfrow=c(3,2), mar=c(3,4,2,2)+.1)

plot(theta, sin2.theta, type='l', ylim=sinRng, xlab='', ylab='Order = 2',
     main='sine(t)' )
lines(theta, sin.theta, lty='dashed')
abline(v=pi3, lty='dotted')

Dsin2.theta = predict(sin2, theta, 1) # first derivative
plot(theta, Dsin2.theta, type='l', ylim=sinRng, xlab='', ylab='',
     main='D sine(t)')
lines(theta, cos(theta), lty='dashed')
abline(v=pi3, lty='dotted')

plot(theta, sin3.theta, type='l', ylim=sinRng, xlab='', ylab='Order = 3')
lines(theta, sin.theta, lty='dashed')
abline(v=pi3, lty='dotted')

Dsin3.theta = predict(sin3, theta, 1)
plot(theta, Dsin3.theta, type='l', ylim=sinRng, xlab='', ylab='')
lines(theta, cos(theta), lty='dashed')
abline(v=pi3, lty='dotted')

plot(theta, sin4.theta, type='l', ylim=sinRng, xlab='t', ylab='Order = 4')
lines(theta, sin.theta, lty='dashed')
abline(v=pi3, lty='dotted')

Dsin4.theta <- predict(sin4, theta, 1)
plot(theta, Dsin4.theta, type='l', ylim=sinRng, xlab='t', ylab='')
lines(theta, cos(theta), lty='dashed')
abline(v=pi3, lty='dotted')

par(op)

#  order 6 spline with equally spaced knots

splinebasis6 = create.bspline.basis(c(0,10), 15, 6)
plot(splinebasis6, lwd=2)

#  order 4 spline with knots not equally spaced

spline.unequal.knots = create.bspline.basis(breaks=c(0, .7, 1))
knots.fda(spline.unequal.knots, interior=FALSE)

#  plot sum of basis function values (all equal to 1)

t10 = seq(0, 10, .1)
spline6.10 = predict(splinebasis6, t10)
plot(t10, rowSums(spline6.10), lwd=2)

#  basis systems for growth data

heightbasis  = with(growth, create.bspline.basis(c(1, 18), 35, 6, age))
plot(heightbasis)

# or

heightbasis. = create.bspline.basis(norder=6, breaks=growth$age)
all.equal(heightbasis, heightbasis.)

help(create.bspline.basis)

##
## Section 3.4 Constant, Monomial and Other Bases
##

conbasis = create.constant.basis(c(0,1))
conb     = create.constant.basis()
all.equal(conbasis, conb)
plot(conbasis)

monbasis = create.monomial.basis(c(0,1), 4)
monb     = create.monomial.basis(nbasis=4)
all.equal(monbasis, monb)
plot(monbasis)

##
## Section 3.5 Methods for Functional Basis Objects
##

methods(class='basisfd')

print(monb)

summary(monb)

monbasis == monb

is(monb, 'basisfd')
inherits(monb, 'basisfd')

monb$params
monb$params = 2*(0:3)

t.1        = seq(0, 1, .1)
monbmat.1  = eval.basis(t.1, monb)
Dmonbmat.1 = eval.basis(t.1, monb, 1)

opar <- par(mfrow=c(1,2))
matplot(monbmat.1, type="l")
matplot(Dmonbmat.1, type="l")
par(opar)

monbmat.1. = predict(monb, t.1)
all.equal(monbmat.1, monbmat.1.)

Dmonbmat.1. = predict(monb, t.1, 1)
all.equal(Dmonbmat.1, Dmonbmat.1.)

##
## Section 3.6 The Structure of the basisfd or basis Class
##

help(basisfd)

