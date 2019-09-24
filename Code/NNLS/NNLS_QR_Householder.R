nnls_QR_House <- function(A, b, ModellpunkteIn = 0, Fehlertol = 1e-10, MainloopMax = 0, InnerloopMax = 1000, StreichloopMax = 1000, Ausgabe = 0, Test = 0) {
  
  N <- dim(A)[2]
  AnzNull <- rowCounts(A, value = 0)
  Zeilen <- (AnzNull != N) & (b != 0) # Wesentliche Zeilen mit Informationen (Zeilen in der nicht alle Werte in A und b gleich 0 sind)
  A <- A[Zeilen,]
  b <- b[Zeilen]
  AnzNull <- AnzNull[Zeilen]
  OrderAnz <- order(AnzNull)
  
  # Initialize
  
  v <- dim(A)[1]
  
  P <- NULL            # At the end, this vector contains the numbers of the columns with positive weighting (= solution)
  n <- 0               # length of P
  R <- 1:N             # At the end, this vector contains the numbers of the columns with weighting zero
  x <- rep(0, N)       # At the end, this vector contains the solution (weighting of all columns)
  Stellen <- -log(Fehlertol, 10)
  
  w <- crossprod(b, A) # calculate gradient "t(A) %*% (b - A %*% x)" for vector x equal to zero 
  wmax <- round(max(w), Stellen)       # maximum gradient
  Abweichung_prev <- round(sum(b^2),5)
  
  Mainloop <- 0
  Streichloop <- 0
  jprev <- 0
  Stabiler_Zustand <- 0
  Abbruch_Grund <- 0
  
  if (ModellpunkteIn == 0) {
    Modellpunkte <- dim(A)[1]
  } else {
    Modellpunkte <- ModellpunkteIn
  }
  Ergebnis_Modellpunkte <- 0
  
  if (MainloopMax == 0) {
    MainloopMax <- dim(A)[2] * 3
  }
  
  ModellpunkteMax <- min (2 * Modellpunkte, dim(A)[1]) # maximum number of model points in the calculation
  
  pb <- winProgressBar(title = "Fortschritt", min = 0, max = ModellpunkteMax, width = 600)
  
  
  # Mainloop
  while ((length(P) < ModellpunkteMax) & (Mainloop < MainloopMax) & (Streichloop < StreichloopMax) & (Stabiler_Zustand == 0) & (wmax > Fehlertol)) {
    Mainloop <- Mainloop + 1
    
    Kandidat_gefunden <- 0
    while (Kandidat_gefunden == 0) {
      j <- R[which.max(w[R])]
      xneu <- A[(n+1):v, j]
      e <- as.matrix(c(1, rep(0, length(xneu)-1)))
      u <- sign(xneu[1]) * sqrt(sum(xneu^2)) * e + xneu
      
      Rlast <- A[(n+1):v,j] - 2 * u %*% (t(u) %*% A[(n+1):v,j]) / as.vector(t(u) %*% u)
      blast <- b[(n+1):v] - 2 * u %*% (t(u) %*% b[(n+1):v]) / as.vector(t(u) %*% u)
      if (Rlast[1] != 0) {
        if (blast[1] / Rlast[1] > 0) {
          Kandidat_gefunden <- 1
        }
      } else {
        print(Mainloop)
        print(j)
        print(Rlast[1])
        print(blast[1])
        print(max(w[R]))
        print(" ")
        w[j] <- 0
      }
    }
    
    j <- which.max(w[R])
    
    if (R[j] == jprev) {
      Stabiler_Zustand <- 1
    }
    
    if (Stabiler_Zustand == 0) {
      
      jprev <- R[j]
      
      if (Ausgabe == 1) {
        print(paste0("Modellpunkt: ",jprev))
      }  
      
      if (length(P) > 0) {
        Rold <- (1:N)[-P]
      } else {
        Rold <- (1:N)
      }
      P <- c(P, jprev)
      
      n <- length(P)
      xneu <- A[n:v, jprev]
      e <- as.matrix(c(1, rep(0, length(xneu)-1)))
      u <- sign(xneu[1]) * sqrt(sum(xneu^2)) * e + xneu
      
      A[n:v,Rold] <- A[n:v,Rold] - 2 * u %*% (t(u) %*% A[n:v,Rold]) / as.vector(t(u) %*% u)
      b[n:v] <- b[n:v] - 2 * u %*% (t(u) %*% b[n:v]) / as.vector(t(u) %*% u)
      
      s <- rep(0, N)
      sP <- matrix(backsolve(A[1:n,P],b), ncol=1)
      s[P] <- sP
      Innerloop <- 0
      
      # Innerloop
      while ((min(sP) <= 0) & (Innerloop < InnerloopMax)) {
        Innerloop <- Innerloop + 1
        neg <- P[sP <= 0]
        alpha <- min (x[neg] / (x[neg] - s[neg]), na.rm = TRUE)
        if ((alpha > 1) | (alpha < 0)) {
          print("Fehler bei alpha-Wert")
        } 
        x[P] <- x[P] + alpha * (s[P] - x[P])
        
        temp <- (round(x[P], Stellen) != 0)
        Spalten <- (1:length(temp))[!temp]
        AnzSpalten <- length(Spalten)
        if (AnzSpalten == 0) {
          print("Fehler => keine Spalte mit 0")
        }
        if (AnzSpalten > 1) {
          print("Mehr als eine Spalte mit 0")
          print(Spalten)
        }
        
        for (LaufSpalten in AnzSpalten:1) {
          if (AnzSpalten == 1) {
            Spalte <- Spalten
          } else {
            Spalte <- Spalten[LaufSpalten]
          }  
          
          P <- P[-Spalte]
          for (Lauf in Spalte:(n-1)) {
            r <- sign(A[Lauf, P[Lauf]]) * sqrt(A[Lauf, P[Lauf]]^2 + A[Lauf+1,P[Lauf]]^2)
            cc <- A[Lauf, P[Lauf]] / r
            ss <- A[Lauf+1, P[Lauf]] / r
            G <- matrix(c(cc,-ss,ss,cc),nrow = 2)
            
            temp <- G %*% A[c(Lauf,Lauf+1),]
            A[c(Lauf,Lauf+1),] <- temp
            
            temp <- G %*% b[c(Lauf,Lauf+1)]
            b[c(Lauf,Lauf+1)] <- temp
            
          }
          n <- length(P)
          
          #          print(Mainloop)
          #          print(Innerloop)
          #          print(Spalte)
          #          print(P)
          
        } # Ende for "LaufSpalten"-Abfrage
        
        sP <- matrix(backsolve(A[1:n,P],b), ncol=1)
        s <- rep(0, N)
        s[P] <- sP
        #        print(s[P])
        #        print(" ")
      } # End Innerloop
      
      Abweichung <- round(sum((b - A[,P] %*% sP)^2),5)
      Max_w_P <- max(crossprod(b - A[,P] %*% sP, A[,P])) * 0
      #      if ((Abweichung > Abweichung_prev) | (Max_w_P > 1)) {
      #        P <- PSicher
      #        RP <- RSicher
      #        QP <- QSicher
      #        QTb <- QTbSicher
      #        x <- xSicher
      #        w <- wSicher
      #        w[jprev] <- 0
      #        Abweichung <- Abweichung_prev
      #        Streichloop <- Streichloop + 1
      #      } else {
      Abweichung_prev <- Abweichung
      x <- s
      w <- crossprod(b - A [,P] %*% sP, A)
      
      if (length(P) == ModellpunkteIn) {
        Ergebnis_Modellpunkte <- 1
        xErg <- x
        PErg <- P
        AErg <- A
        bErg <- b
        wErg <- w
      }
      
      R <- (1:N)[-P]
      wmax <- round(max(w[R]), Stellen)
      
      if (Test == 0) {      
        setWinProgressBar(pb, length(P), title = paste0( "Modellpunkte ", length(P), " von max. ", Modellpunkte, " Mainloop: ", Mainloop))
      } else if (Test == 1) {
        setWinProgressBar(pb, length(P), title = paste0( "Modellpunkte ", length(P), " von max. ", Modellpunkte, " Mainloop: ", Mainloop, " Diff: ", round(Abweichung, 0), " max(w): ", round(wmax, 0), " Kand: ", sum(w[R] > Fehlertol)))
      }  
    } # End Stabiler_Zustand
  } # End Mainloop
  
  close(pb)
  
  if (Ergebnis_Modellpunkte == 1) {
    P <- PErg
    x <- xErg
    A <- AErg
    b <- bErg
    w <- wErg
    R <- (1:N)[-P]
  }
  
  fitted <- A %*% x
  resid <-  b - fitted
  if (Stabiler_Zustand == 1) {
    Abbruch_Grund <- 1
    print("Stabiler Zustand erreicht")
  }
  else if (Mainloop >= MainloopMax) {
    print("Maximale Anzahl der Mainloop erreicht")
    Abbruch_Grund <- 2
  }
  else if (Streichloop >= StreichloopMax) {
    print("Maximale Anzahl der Streichloop erreicht")
    Abbruch_Grund <- 3
  }
  else if (max(w) <= Fehlertol) {Abbruch_Grund <- 0}
  
  xepsilon <- x
  Pepsilon <- (1:N)[xepsilon > Fehlertol]
  xepsilon[xepsilon <= Fehlertol] <- 0
  
  nnls_QR.out <- list( x = x,
                       A = A,
                       b = b,
                       xepsilon = xepsilon,
                       deviance = sum(resid^2),
                       residuals = resid, 
                       fitted = fitted,
                       mw = max(w),
                       w = w,
                       passiveorig = P,
                       passive = sort(P),
                       Pepsilon = Pepsilon,
                       loops = Mainloop,
                       nsetp = length(P),
                       Abbruch = Abbruch_Grund)
  
  class(nnls_QR.out) <- "nnls"
  nnls_QR.out 
} # End nnls_QR_House




