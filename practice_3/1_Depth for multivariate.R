#https://cran.r-project.org/web/packages/DepthProc/DepthProc.pdf
#http://www.csun.edu/~ctoth/Handbook/chap58.pdf
library(DepthProc)
library(MultiRNG)

### A simple example
x <- c(1, 3, 7, 9, 15)
x
median(x)

dp <- matrix(NA, nrow = 5, ncol = 4)
rownames(dp) <- as.factor(x)
colnames(dp) <- c("Projection", "Mahalanobis", "Euclidean", "Tukey")

dp[,1] <- depth(x, x, method = "Projection")
dp[,2] <- depth(x, x, method = "Mahalanobis")
dp[,3] <- depth(x, x, method = "Euclidean")
dp[,4] <- depth(x, x, method = "Tukey")

dp


y <- c(1, 3, 7, 9, 15, 19)
y
median(y)

dp2 <- matrix(NA, nrow = 6, ncol = 4)
rownames(dp2) <- as.factor(y)
colnames(dp2) <- c("Projection", "Mahalanobis", "Euclidean", "Tukey")

dp2[,1] <- depth(y, y, method = "Projection")
dp2[,2] <- depth(y, y, method = "Mahalanobis")
dp2[,3] <- depth(y, y, method = "Euclidean")
dp2[,4] <- depth(y, y, method = "Tukey")

dp2 # 7 and 9 are the deepest points but not for all the depths.
# ideed accordi to Euclidean the point 9 is the deepest 

### Univariate case

# Compute depth of a point with respect to a data set
set.seed(123)
w1 <- runif(200, 0, 1)
depth(mean(w1), w1)
depth(median(w1), w1)
median(w1)
sort(depth(w1,w1),decreasing=TRUE) 

# Find maximal depth
w2 <- rnorm(200, 0, 1)
depth_val <- numeric(200)
for(i in 1:200) {
  depth_val[i] <- depth(w2[i], w2)
}
#maximal depth
max(depth_val)
#which
w2[which.max(depth_val)]
median(w2)

####
data(starsCYG, package = "robustbase")
head(starsCYG)

st <- as.matrix(starsCYG)
dim(st)

plot(starsCYG)

# Euclidean depth
dE <- depthEuclid(st, st)
mdE <- which.max(dE)
dE[mdE]
st[mdE,] # deepest point

apply(st,2,median)

plot(st, xlim=c(3,5), ylim=c(3,7))
points(st[mdE,1], st[mdE,2], pch=23, col="blue", bg="blue", lwd=2)
points(median(st[,1]), median(st[,2]), pch=24, col="red", bg="red", lwd=2)

# red one is median point while the blue one is the deepest point, affected by outliers

# Local depth
dL <- depthLocal(st, st, depth_params1 = list(method = "LP"))
dL
mdL <- which.max(dL)
dL[mdL]
st[mdL,]

depthContour(st, depth_params = list(method = "Local", depth_params1 = list(method = "LP"))) # 3D plot

# MBD, Frainman-Muniz
dMBD <- fncDepth(st, method = "MBD")
dFM <- fncDepth(st, method = "FM")
mdMBD <- which.max(dMBD)
mdFM <- which.max(dFM)

dMBD[mdMBD]
st[mdMBD,]

dFM[mdFM]
st[mdFM,]

plot(st, xlim=c(3,5), ylim=c(3,7))
points(st[mdMBD,1], st[mdMBD,2], pch=23, col="blue", bg="blue", lwd=2)
points(st[mdFM,1], st[mdFM,2], pch=23, col="green", bg="green", lwd=2)
points(median(st[,1]), median(st[,2]), pch=24, col="red", bg="red", lwd=2)

fncDepthMedian(st, method = "MBD")
fncDepthMedian(st, method = "FM")

