## ----echo=FALSE,include=FALSE-------------------------------------------------
library(QR)
set.seed(26112011)

## ----echo=TRUE, include=TRUE--------------------------------------------------
# Let's sample a random square-matrix
A<-matrix(runif(121,min = -100, max = 100), 11, 11)

## ----echo=TRUE, include=TRUE, collapse=TRUE-----------------------------------
QRres<-QR(A)
QRres

## ----echo=TRUE, include=TRUE, collapse=TRUE-----------------------------------
QRres$Q

QRres$R

## ----echo=TRUE, include=TRUE, collapse=TRUE-----------------------------------
all.equal(QRres$Q%*%t(QRres$Q),diag(11))

## ----echo=TRUE, include=TRUE, collapse=TRUE-----------------------------------
all.equal(QRres$R[lower.tri(QRres$R)],rep(0,11*10/2))

## ----echo=TRUE, include=TRUE, collapse=TRUE-----------------------------------
all.equal(QRres$Q%*%QRres$R,A)

## ----echo=TRUE, include=TRUE, collapse=TRUE-----------------------------------
# example of pivoting
x <- cbind(int = 1,
           b1 = rep(1:0, each = 3), b2 = rep(0:1, each = 3),
           c1 = rep(c(1,0,0), 2), c2 = rep(c(0,1,0), 2), c3 = rep(c(0,0,1),2))
x # is singular, columns "b2" and "c3" are "extra"
a <- qr(x)

## ----echo=TRUE, include=TRUE, collapse=TRUE-----------------------------------
all.equal(qr.Q(a)%*%qr.R(a),x)

## ----echo=TRUE, include=TRUE, collapse=TRUE-----------------------------------
all.equal(QR(x)$Q%*%QR(x)$R,x)

