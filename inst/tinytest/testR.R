#Check R is upper triangular
#Low triangle is 0.

#Square matrix
A<-matrix(sample(-100:100, size = 49, replace = TRUE),7,7)
R<-QR(A)$R

expect_equal(R[lower.tri(R)],rep(0,6*7/2))



#(m x n) matrix with m > n.

B<-matrix(sample(-100:100, size = 35, replace = TRUE),7,5)
Rb<-QR(B)$R


expect_equal(Rb[lower.tri(Rb)],rep(0,4*5/2))



#(m x n) matrix with m < n

C<-matrix(sample(-100:100, size = 35, replace = TRUE),5,7)
Rc<-QR(C)$R

expect_equal(Rc[lower.tri(Rc)],rep(0,4*5/2))
