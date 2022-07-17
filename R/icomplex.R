#' QR factorization without pivoting
#'
#' This function performs QR factorization without pivoting to a complex matrix A.
#'
#'
#' @param A a complex matrix whose QR decomposition is to be computed.
#' @param complete boolean that indicates if the R matrix should be completed with 0s
#' to its full rank.
#'
#' @details This method is an alternative to the default qr function of base R.
#' The default function returns a pivoted solution in many cases, which is
#' not always the desired solution. In this function, we returned the unpivoted
#' solution for the QR factorization using the LAPACK routine ZGEQRF. This function
#' works with real and complex matrices.
#'
#' @return Returns a list with the following components:
#'  \item{qr}{a matrix with the same dimensions as A. The upper triangle contains
#'  the \bold{R} of the decomposition and the lower triangle contains information
#'  on the \bold{Q} of the decomposition (stored in compact form).}
#'  \item{qraux}{a vector of length ncol(A) which contains additional
#'  information on \bold{Q}.}
#'  \item{Q}{an orthogonal matrix such that Q*R is the input matrix.}
#'  \item{R}{an upper triangular matrix such that Q*R is the input matrix.}
#'
#'
#' @source
#'
#' LAPACK routine ZGEQRF is used for the QR factorization without pivoting.
#'
#'
#' @references Anderson. E. and ten others (1999) LAPACK Users' Guide. Third Edition. SIAM.
#' Available on-line at http://www.netlib.org/lapack/lug/lapack_lug.html.
#'
#'
#' @examples
#'
#' set.seed(2)
#' A<-matrix(c(complex(real=1,imaginary = 1),
#' complex(real=3,imaginary = -2), complex(real=2,imaginary = 1),
#' complex(real=0,imaginary = 3)),2,2)
#' qres<-QRcomp(A)
#'
#' #Inspect the main results of the factorization:
#' qres$Q
#' qres$R
#'
#' @keywords QR
#' @keywords factorization
#'
#' @export
#'
#' @useDynLib QR, .registration=TRUE
QRcomp <- function(A, complete=FALSE)
{

  if(is.numeric(A)){
    nr<-NROW(A)
    nc<-NCOL(A)
    A<-matrix(as.complex(A),nrow = nr, ncol = nc)
  }

  stopifnot(is.complex(A))
  stopifnot(is.finite(A))
  stopifnot(is.matrix(A))

  qrC<-.C("qrcomplexc", MX = nrow(A), NX = ncol(A), X=A, # X=matrix(as.complex(A),nrow(A),ncol(A)),
          LDX=nrow(A), tau=rep(as.complex(0),min(nrow(A),ncol(A))), outlwork=as.integer(1), PACKAGE = "QR")

  qrC$qr<-qrC$X

  R<-qrRcomp(qrC, complete)

  Q<-qrQcomp(qrC$X,qrC$tau)

  out<-list(qr=qrC$X, qraux=qrC$tau, Q=Q, R=R)
  return(out)
}



qrRcomp<-function(qr, complete = FALSE){
  R <- qr$X
  if(!complete){
    R <- R[seq.int(min(dim(R))), , drop = FALSE]
  }
  R[row(R) > col(R)] <- 0
  return(R)
}


qrQcomp<-function(qr,tau){
  k<-min(NROW(qr),NCOL(qr))
  #Id<-diag(nrow = k)
  k2<-max(NROW(qr),NCOL(qr))

  if(NROW(qr)<NCOL(qr)){
    Id<-diag(nrow = k)
  }else{
    Id<-diag(nrow = k2)
  }


  Q<-Id
  for (idx in 1:k) {
    if(idx==1){
      v<-matrix(c(1,qr[2:NROW(qr),idx]),ncol = 1)
    }else if(idx==NROW(qr) && NROW(qr)<=NCOL(qr)){
      v<-matrix(c(rep(0,idx-1),1),ncol = 1)
    }else{
      v<-matrix(c(rep(0,idx-1),1,qr[(idx+1):NROW(qr),idx]),ncol = 1)
    }
    H<-Id-tau[idx]*v%*%t(Conj(v))
    Q<-Q%*%H
  }

  if(NROW(qr)<NCOL(qr)){
    Q<-Q[1:NROW(qr),1:NROW(qr)]
  }else{
    Q<-Q[1:NROW(qr),1:NCOL(qr)]
  }
  return(Q)
}



