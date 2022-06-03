#Check that Q*R is equal to the original matrix

#Square matrix
A<-matrix(sample(-100:100, size = 49, replace = TRUE),7,7)

expect_equal(QR(A)$Q%*%QR(A)$R,A)



#(m x n) matrix with m > n.

B<-matrix(sample(-100:100, size = 35, replace = TRUE),7,5)

expect_equal(QR(B)$Q%*%QR(B)$R,B)



#(m x n) matrix with m < n

C<-matrix(sample(-100:100, size = 35, replace = TRUE),5,7)

expect_equal(QR(C)$Q%*%QR(C)$R,C)
