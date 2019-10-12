# The pseudocode for the nnls-algoritm can be found at
# https://en.wikipedia.org/wiki/Non-negative_least_squares
nnls_VMF <- function(A, b, Modellpunkte = 0, Fehlertol = 1e-10, MainloopMax = 10000, InnerloopMax = 1000, Ausgabe = 0) {
  
  N <- dim(A)[2]
  # Count the number of zeros per row; if for one row there are N zeros all policies
  # produce zeros for that cashflow --> not a relevant information. 
  AnzNull <- rowCounts(A, value = 0)
  # find those rows which have non-zero values in A and b simultaneously
  # Those rows have relevant information
  Zeilen <- (AnzNull != N) & (b != 0)
  # Reduce A and b to those rows with relevant information
  A <- A[Zeilen,]
  b <- b[Zeilen]
  
  # Initialize
  P <- NULL
  R <- 1:N
  x <- rep(0, N)
  
  # Log the number of non-zero entries in x by counting the elements in P
  nonZeroEntries <- length(P)
  
  # Log the development of the residuals
  fitted <- A %*% x
  residual_tmp <-  (b - fitted)^2 %>% sum()
  
  # Calculate initial Gradient "t(A) %*% (b - A %*% x)" 
  # Since in first iteration x is zero this simplifies to "t(A) %*% b"
  w <- crossprod(b, A)  
  # Find highest descent 
  wmax <- max(w)
  Stellen <- -log(Fehlertol, 10)
  # initialize counter for Main loop
  Mainloop <- 0
  
  # initialise variables for the exit-reason
  jprev <- 0
  Stabiler_Zustand <- 0
  Abbruch_Grund <- 0
  
  # define the tolerance based on the machine precision 
  TestTol <- .Machine$double.eps * 1000
  
  # When the user doesn't restrict the number of modelpoints then 
  # we set the number of maximal modelpoints to the number of varibles we 
  # want to match. 
  if (Modellpunkte == 0) {
    Modellpunkte <- dim(A)[1]
  }
  ModellpunkteMax <- min(2*Modellpunkte, dim(A)[1])
  
  # initialize progress bar
  pb <- winProgressBar(title = "Progress",
                       min = 0,
                       max = Modellpunkte,
                       width = 600)
  
  # Mainloop
  while ((length(P) < ModellpunkteMax) &
         (Mainloop < MainloopMax) & 
         (Stabiler_Zustand == 0) &
         (wmax > Fehlertol)) {
    
    # log development of residuals
    fitted <- A %*% x
    resid <-  b - fitted

    
    Stabiltest <- 0
    while (Stabiltest == 0) {
      j <- which.max(w[R])
      Kandidat <- R[j] 
      Ptemp <- unique(sort(rbind(P, Kandidat)))
      # check if det below tolerance, if this is the case then we can't invert
      # t(A) %*% A restricted to the columns of P
      if (det(crossprod(A[,Ptemp])) <= TestTol) {
        w[Kandidat] <- 0
      } else {
        Stabiltest <- 1
      }  
    }
    
    # if the new candidate is the same as in the last iteration of the loop
    # a stable state is reached and the optimization is finished
    if (R[j] == jprev) {
      #Stabiler_Zustand <- 1
    }
    
    #if (Stabiler_Zustand == 0) {
    if (0 == 0) {
      jprev <- R[j]
      
      # move index form R to P
      P <- unique(sort(rbind(P, R[j])))
      
      # define restricted matrix A^P
      AP <- A[, P]
      s <- rep(0, N)
      # compute the restricted least squres problem
      sP <- solve( crossprod(AP), crossprod(AP, b))
      #      sP <- solve( t(AP) %*% AP, t(AP) %*% b)
      #      sP <- as.matrix(lsfit( AP, b, intercept = FALSE)[[1]])
      # assign the results of the restriced least squres problem to the 
      # entries of s indexed in P.
      s[P] <- sP
      Innerloop <- 0
      
      # Innerloop
      while ((min(sP) <= 0) & (Innerloop < InnerloopMax)) {
        # find indices for negative values
        neg <- P[sP <= 0]
        alpha <- min (x[neg] / (x[neg] - s[neg]), na.rm = TRUE)
        x[P] <- x[P] + alpha * (s[P] - x[P])

        temp <- (round(x[P], Stellen) != 0)
        AP <- AP[,temp]
        P <- P[temp]
        sP <- solve( crossprod(AP), crossprod(AP, b))
        #        sP <- solve( t(AP) %*% AP, t(AP) %*% b)
        #        sP <- as.matrix(lsfit( AP, b, intercept = FALSE)[[1]])
        s <- rep(0, N)
        s[P] <- sP
        Innerloop <- Innerloop + 1
      } # End Innerloop
 
      R <- (1:N)[-P]
      x <- s
      w <- crossprod(b - AP %*% sP, A)
      wmax <- max(w)
      Mainloop <- Mainloop + 1
      
      
      # logging
      nonZeroEntries <- rbind(nonZeroEntries, length(P))
      # Log the development of the residuals
      fitted <- A %*% x
      residual_tmp <-  c(residual_tmp, (b - fitted)^2 %>% sum())
      
      
      
      setWinProgressBar(
        pb,
        length(P),
        title = paste0("Policies ", length(P), " of max. ",Modellpunkte,
                       " Mainloop: ", Mainloop, " max(w): ", round(wmax, 4)
        )
      )
    } # End Stabiler_Zustand
  } # End Mainloop
  
  # close progress bar
  close(pb)
  
  # Apply fit and calculate residuals
  fitted <- A %*% x
  resid <-  b - fitted
  
  if (Stabiler_Zustand == 1) {
    Abbruch_Grund <- "stable state"
  }
  else if (Mainloop >= MainloopMax) {
    Abbruch_Grund <- "main loop max reached"
  }
  else if (max(w) <= Fehlertol) {
    Abbruch_Grund <- "gradient too small"
  }
  
  xepsilon <- x
  Pepsilon <- (1:N)[xepsilon > Fehlertol]
  xepsilon[xepsilon <= Fehlertol] <- 0
  
  nnls_VMF.out <- list(
    x = x,
    xepsilon = xepsilon,
    deviance = sum(resid ^ 2),
    residuals = resid,
    fitted = fitted,
    mw = max(w),
    w = w,
    passive = P,
    Pepsilon = Pepsilon,
    loops = Mainloop,
    nsetp = length(P),
    Abbruch = Abbruch_Grund,
    Entwicklung = nonZeroEntries,
    Residual_Entwicklung = residual_tmp
  )
  
  class(nnls_VMF.out) <- "nnls"
  nnls_VMF.out 
} # End nnls_VMF