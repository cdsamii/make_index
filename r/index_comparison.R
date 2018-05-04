rm(list=ls())

###################################################
# Functions:

# Function to standardize columns of a matrix
# where you designate a standardization group
# (e.g., the control group in an experiment)
# with "sgroup", a logical vector.

matStand <- function(x, sgroup = rep(TRUE, nrow(x))){
				for(j in 1:ncol(x)){
					x[,j] <- (x[,j] - mean(x[sgroup,j]))/sd(x[sgroup,j])
				}
				return(x)
			}

# Function that takes in data in matrix format and returns
# (i) IC weights and (ii) ICW index.
# Weights can be incorporated using the "wgts" argument.
# The "revcols" argument takes a vector indicating which columns,
# if any, should have values reversed (that is, standardized 
# values are multiplied by -1) prior to construction of the index. 

icwIndex <- function(	xmat,
						wgts=rep(1, nrow(xmat)),
						revcols = NULL,
						sgroup = rep(TRUE, nrow(xmat))){
					X <- matStand(xmat, sgroup)
					if(length(revcols)>0){
						X[,revcols] <-  -1*X[,revcols]
					}
					i.vec <- as.matrix(rep(1,ncol(xmat)))
					Sx <- cov.wt(X, wt=wgts)[[1]]
					weights <- solve(t(i.vec)%*%solve(Sx)%*%i.vec)%*%t(i.vec)%*%solve(Sx)
					index <- t(solve(t(i.vec)%*%solve(Sx)%*%i.vec)%*%t(i.vec)%*%solve(Sx)%*%t(X))
					return(list(weights = weights, index = index))
					}

# Factor analysis function (R base does not have one. factanal() is MLE, which we don't want.) 
source("http://www.stat.sc.edu/~habing/courses/530/fact.txt")
###################################################

# Inverse covariance weighted index

# Nothing is correlated case
n <- 2000

x1.l <- rnorm(n)
x1.raw <- as.numeric(cut(x1.l, c((min(x1.l)-1),-1,0,1,(max(x1.l)+1)), labels=c(1,2,3,4)))
x2.raw <- rbinom(n,1,.25)
x3.l <- rnorm(n)
x3.raw <- as.numeric(cut(x3.l, c((min(x3.l)-1),-1,0,1,(max(x3.l)+1)), labels=c(1,2,3,4)))

X.raw <- cbind(x1.raw,x2.raw,x3.raw)

icwX <- icwIndex(X.raw)

# The weights that are applied:
cor(X.raw)
icwX$weights

# The resulting index

z <- icwX$index

mydata <- data.frame(z,jitter(x1.raw),jitter(x2.raw),jitter(x3.raw))

par(mar=c(0,0,0,0)+4)
pairs(mydata)

###################################################

# Now, make x1 and x3 correlated (x2 remains orthogonal):
x1 <- (x1.raw - mean(x1.raw))/sd(x1.raw)
x3.l <- rnorm(n,mean=2*x1)
x3.raw <- as.numeric(cut(x3.l, c((min(x3.l)-1),-1,0,1,(max(x3.l)+1)), labels=c(1,2,3,4)))
X <- cbind(x1.raw,x2.raw,x3.raw)
cor(X)

# Inverse covariance weighting

icwX.2 <- icwIndex(X)
icwX.2$weights
z <- icwX.2$index

mydata <- data.frame(z,jitter(x1.raw),jitter(x2.raw),jitter(x3.raw))
pairs(mydata)

# Now compare to factor analysis and principal components

# Factor Analysis

X.raw <- as.matrix(data.frame(x1.raw,x2.raw,x3.raw))

lambda <- fact(X.raw,method="iter",maxfactors=2, niter=130)$loadings
sigma <- cor(X.raw)
# This Thompson's method.
z.fac <- matStand(X.raw)%*%solve(sigma)%*%lambda

# Compare the weightings: 

# ICW: should see ca. 25% for x1 and x3, 50% for x2:
icwX.2$weights

# Factor analysis: factor 1 loads x1 and x2 equally, no weight to x2; 
# vice versa for second factor: 
t(lambda)

z.mat <- data.frame(z.icw = z, z.fac = z.fac, x1.raw,x2.raw,x3.raw)
pairs(z.mat)