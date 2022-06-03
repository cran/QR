#Check Q is orthogonal
#Q%*%t(Q)==identityMatrix

#Square matrix
A<-matrix(sample(-100:100, size = 49, replace = TRUE),7,7)
Q<-QR(A)$Q
Id<-diag(nrow = 7, ncol = 7)

expect_equal(Q%*%t(Q),Id)



#(m x n) matrix with m > n. Issue here. Q same as in Matlab and qr(). In matlab orthogonal and here not.

B<-matrix(sample(-100:100, size = 35, replace = TRUE),7,5)
Qb<-QR(B)$Q
Idb<-diag(nrow = 5, ncol = 5)

expect_equal(t(Qb)%*%Qb,Idb)




#(m x n) matrix with m < n

C<-matrix(sample(-100:100, size = 35, replace = TRUE),5,7)
Qc<-QR(C)$Q
Idc<-diag(nrow = 5, ncol = 5)

expect_equal(Qc%*%t(Qc),Idc)
