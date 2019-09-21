#' ---
#' Title:   "The file contains some compression algorithm (Cluster, Optimisation)"
#' Author:  ""
#' Date:    "July 2019"
#' ---

# The pseudocode for the nnls-algoritm can be found in https://en.wikipedia.org/wiki/Non-negative_least_squares
nnls_VMF <- function(A, b, Modellpunkte = 0, Fehlertol = 1e-10, MainloopMax = 10000, InnerloopMax = 1000, Ausgabe = 0) {
  
  N <- dim(A)[2]
  AnzNull <- rowCounts(A, value = 0)
  Zeilen <- (AnzNull != N) & (b != 0) # Wesentliche Zeilen mit Informationen (Zeilen in der nicht alle Werte in A und b gleich 0 sind)
  A <- A[Zeilen,]
  b <- b[Zeilen]
  AnzNull <- AnzNull[Zeilen]
  
  # Initialize
  
  N <- dim(A)[2]
  P <- NULL
  R <- 1:N
  x <- rep(0, N)
  
  nonZeroEntries <- c(0)
  
  # calculate initial Gradient "t(A) %*% (b - A %*% x)" 
  # since in first iteration x is zero this simplifies to "t(A) %*% b"
  w <- crossprod(b, A)  
  # find highest descent 
  wmax <- max(w)
  Stellen <- -log(Fehlertol, 10)
  # initialize counter for Main loop
  Mainloop <- 0
  
  #  InnerloopSum <- 0
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
  
  
  pb <- winProgressBar(title = "Fortschritt", min = 0, max = Modellpunkte, width = 600)
  
  # Mainloop
  while ((length(P) < ModellpunkteMax) & (Mainloop < MainloopMax) & (Stabiler_Zustand == 0) & (wmax > Fehlertol)) {
    
    # 
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
    
    if (R[j] == jprev) {
      Stabiler_Zustand <- 1
    }
    
    if (Stabiler_Zustand == 0) {
      jprev <- R[j]
      #      print(jprev)
      
      if (Ausgabe == 1) {
        print(paste0("Modellpunkt: ",R[j]))
      }  
      P <- unique(sort(rbind(P, R[j])))
      #      R <- (1:N)[-P]
      
      AP <- A[,P]
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
        neg <- P[sP <= 0]
        alpha <- min (x[neg] / (x[neg] - s[neg]), na.rm = TRUE)
        x[P] <- x[P] + alpha * (s[P] - x[P])
        #        temp <- which.min(x[P])
        #        AP <- AP[,-temp]
        #        P <- P[-temp]
        temp <- (round(x[P], Stellen) != 0)
        AP <- AP[,temp]
        P <- P[temp]
        sP <- solve( crossprod(AP), crossprod(AP, b))
        #        sP <- solve( t(AP) %*% AP, t(AP) %*% b)
        #        sP <- as.matrix(lsfit( AP, b, intercept = FALSE)[[1]])
        s <- rep(0, N)
        s[P] <- sP
        Innerloop <- Innerloop + 1
        #        InnerloopSum <- InnerloopSum + 1
        #        setWinProgressBar(pb, length(P), title = paste0( "Modellpunkte ", length(P), " von max. ", Modellpunkte, " Mainloop: ", Mainloop, " Innerloop: ", Innerloop, " max(w): ", round(wmax, 4)))
        
      } # End Innerloop
      if (Ausgabe == 1) {
        print(paste0("Anzahl innere Schleifen: ",Innerloop))
      }  
      R <- (1:N)[-P]
      x <- s
      w <- crossprod(b - AP %*% sP, A)
      wmax <- max(w)
      Mainloop <- Mainloop + 1
      
      nonZeroEntries <- rbind(nonZeroEntries, length(P))
      
      setWinProgressBar(pb, length(P), title = paste0( "Modellpunkte ", length(P), " von max. ", Modellpunkte, " Mainloop: ", Mainloop, " max(w): ", round(wmax, 4)))
    } # End Stabiler_Zustand
  } # End Mainloop
  
  close(pb)
  
  fitted <- A %*% x
  resid <-  b - fitted
  if (Stabiler_Zustand == 1) {Abbruch_Grund <- 1}
  else if (Mainloop >= MainloopMax) {
    print("Maximale Anzahl der Mainloop erreicht")
    Abbruch_Grund <- 2
  }
  else if (max(w) <= Fehlertol) {Abbruch_Grund <- 3}
  xepsilon <- x
  Pepsilon <- (1:N)[xepsilon > Fehlertol]
  xepsilon[xepsilon <= Fehlertol] <- 0
  
  nnls_VMF.out <- list( x = x,
                        xepsilon = xepsilon,
                        deviance = sum(resid^2),
                        residuals = resid, 
                        fitted = fitted,
                        mw = max(w),
                        w = w,
                        passive = P,
                        Pepsilon = Pepsilon,
                        loops = Mainloop,
                        nsetp = length(P),
                        Abbruch = Abbruch_Grund,
                        Entwicklung = nonZeroEntries)
  
  class(nnls_VMF.out) <- "nnls"
  nnls_VMF.out 
} # End nnls_VMF