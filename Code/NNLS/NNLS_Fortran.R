nnls_Fortran <- function(A, b) {
  MDA <- M <- nrow(A) 
  N <- ncol(A)
  RNORM <- MODE <- 0
  NSETP <- 0
  W <- INDEX <- X <- rep(0, N)
  ZZ <- rep(0, M)
  sol <- .Fortran("nnls", A = as.numeric(A), MDA = as.integer(MDA), M =
                    as.integer(M), N = as.integer(N), B = as.numeric(b),
                    X = as.numeric(X), RNORM = as.numeric(RNORM), W =
                    as.numeric(W), ZZ = as.numeric(ZZ), INDEX =
                    as.integer(INDEX), MODE = as.integer(MODE),
                    NSETP = as.integer(NSETP), PACKAGE="nnls")
  fitted <- A %*% sol$X
  resid <-  b - fitted
  index <- sol$INDEX
  nsetp <- sol$NSETP
  if(nsetp > 0)
    passive <- index[1:nsetp]
  else passive <- vector()
  if(nsetp == N)
    bound <- vector() 
  else 
    bound <- index[(nsetp+1):N]
  nnls.out <- list(x=sol$X,
                   A=sol$A,
                   b=sol$B,
                   w=sol$W,
                   deviance=sol$RNORM^2,
                   residuals=resid, 
                   fitted = fitted,
                   mode=sol$MODE,
                   passive = passive, 
                   bound = bound, 
                   nsetp = nsetp)
  class(nnls.out) <- "nnls"
  nnls.out 
} # End nnls_Fortran
