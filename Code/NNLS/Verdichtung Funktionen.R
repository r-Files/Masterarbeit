#' ---
#' Title:   "The file contains some compression algorithm (Cluster, Optimisation)"
#' Author:  ""
#' Date:    "July 2019"
#' ---

# The pseudocode for the nnls-algoritm can be found in https://en.wikipedia.org/wiki/Non-negative_least_squares

nnls_QR <- function(A, b, ModellpunkteIn = 0, Fehlertol = 1e-10, MainloopMax = 0, InnerloopMax = 1000, StreichloopMax = 1000, 
                    Ausgabe = 0, Test = 0, Uebergabe = 0, PIn, RPIn, QPIn, QTbIn, xin) {
  
  N <- dim(A)[2]
  AnzNull <- rowCounts(A, value = 0)
  Zeilen <- (AnzNull != N) & (b != 0) # Wesentliche Zeilen mit Informationen (Zeilen in der nicht alle Werte in A und b gleich 0 sind)
  A <- A[Zeilen,]
  b <- b[Zeilen]
  AnzNull <- AnzNull[Zeilen]

  # Initialize

  N <- dim(A)[2]
  P <- NULL            # At the end, this vector contains the numbers of the columns with positive weighting (= solution)
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

  if (Uebergabe == 1) {
    P <- PIn
    R <- (1:N)[-P]
    RP <- RPIn
    QP <- QPIn
    QTb <- QTbIn
    x[P] <- xin[P]
    w <- crossprod(b - A[,P] %*% x[P], A)
    wmax <- round(max(w[R]), Stellen)
    Mainloop <- 1
    print("Uebergabe")
    print(length(P))
  }
    
  # Mainloop
  while ((length(P) < ModellpunkteMax) & (Mainloop < MainloopMax) & (Streichloop < StreichloopMax) & (Stabiler_Zustand == 0) & (wmax > Fehlertol)) {
    Mainloop <- Mainloop + 1

    # Sicherung aktueller Daten
    if ((Mainloop > 1) & (length(P) > 1)) {
      PSicher <- P
      RSicher <- RP
      QSicher <- QP
      QTbSicher <- QTb
      xSicher <- x
      wSicher <- w
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
      P <- c(P, jprev)
      AP <- A[,P]
      
      if (length(P) == 1){ # decomposition of a matrix with only one column
        RP <- as.double(sqrt(t(AP) %*% AP))
        QP <- AP / RP
        QTb <- t(QP) %*% b
      } else {
        u <- A[,jprev]
        u <- u - QP %*% t(u %*% QP)
        unorm <- as.double(sqrt(t(u) %*% u))
        u <- u / unorm
        RP <- cbind(rbind(RP,0), c(A[,jprev] %*% QP, unorm))
        QP <- cbind(QP, u)
        QTb <- c(QTb, t(u) %*% b)
        if (sign(RP[length(P),length(P)]) != sign(QTb[length(P)])) {
          print("Unterschiedliches VZ")
        }
      }
      
      s <- rep(0, N)
      sP <- matrix(backsolve(RP, QTb), ncol=1)
      s[P] <- sP
      Innerloop <- 0
      
      # Innerloop
#      xsicher <- x
      while ((min(sP) <= 0) & (Innerloop < InnerloopMax)) {
        Innerloop <- Innerloop + 1
#        x <- xsicher
        neg <- P[sP <= 0]
        alpha <- min (x[neg] / (x[neg] - s[neg]), na.rm = TRUE)
        if ((alpha >= 1) | (alpha <= 0)) {
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
          if (1 == 2) { # Soll die Givens-Funktion als Funktion verwendet werden => wird dadurch langsamer
            ErgGiv <- Givens_VMF(QP, RP, QTb, Spalte)

            QP <- ErgGiv$Q
            RP <- ErgGiv$R
            QTb <- ErgGiv$QTb
          } else {
          
            if (Spalte == (dim(RP)[2])) { # Die letzte Spalte wird gestrichen => nichts zu berechnen
              RNeu <- RP
              QNeu <- QP
              QTbNeu <- QTb
              print("Letzte Spalte wird gestrichen")
            } else {
              n <- dim(QP)[2]
              RNeu <- RP
              QNeu <- QP
              QTbNeu <- QTb
              
              for (Lauf in (Spalte+1):n) { # Givens-Transformation
                r <- sign(RNeu[Lauf, Lauf]) * sqrt(RNeu[Lauf, Lauf]^2 + RNeu[Spalte,Lauf]^2)
                cc <- RNeu[Lauf, Lauf] / r
                ss <- RNeu[Spalte, Lauf] / r
                G <- matrix(c(cc,-ss,ss,cc),nrow = 2)
                
                v <- Lauf:n
                temp <- G %*% RNeu[c(Lauf,Spalte),v]
                RNeu[c(Lauf,Spalte),v] <- temp
                
                temp <- QNeu[,c(Lauf,Spalte)] %*% t(G)
                QNeu[,c(Lauf,Spalte)] <- temp
                
                temp <- G %*% QTbNeu[c(Lauf,Spalte)]
                QTbNeu[c(Lauf,Spalte)] <- temp
              } # Ende for-Schleife "Lauf
            } # Ende if "Spalte"-Abfrage
            QP <- QNeu[,-Spalte]
            RP <- RNeu[-Spalte,-Spalte]
            QTb <- QTbNeu[-Spalte]
          } # Ende if Givens
          AP <- AP[,-Spalte]
          P <- P[-Spalte]
        } # Ende for "LaufSpalten"-Abfrage
        
        sP <- matrix(backsolve( RP, QTb), ncol=1)
        s <- rep(0, N)
        s[P] <- sP
      } # End Innerloop

      if (Test == 1) {
        if (length(P) > 2) {
          if (max(abs(AP - QP %*% RP)) > 1e-5) {
            print("Fehler")
          }
        }
      } # Ende if "Test"

      Abweichung <- round(sum((b - AP %*% sP)^2),5)
      Max_w_P <- max(crossprod(b - AP %*% sP, AP)) * 0
      if ((Abweichung > Abweichung_prev) | (Max_w_P > 1)) {
        P <- PSicher
        RP <- RSicher
        QP <- QSicher
        QTb <- QTbSicher
        x <- xSicher
        w <- wSicher
        w[jprev] <- 0
        Abweichung <- Abweichung_prev
        Streichloop <- Streichloop + 1
      } else {
        Abweichung_prev <- Abweichung
        x <- s
        w <- crossprod(b - AP %*% sP, A)
        Streichloop <- 0

        if (length(P) == ModellpunkteIn) {
          Ergebnis_Modellpunkte <- 1
          xErg <- x
          PErg <- P
          APErg <- AP
          QPErg <- QP
          RPErg <- RP
          QTbErg <- QTb
          wErg <- w
        }
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
    AP <- APErg
    QP <- QPErg
    RP <- RPErg
    QTb <- QTbErg
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
                       AP = AP,
                       b = b,
                       QP = QP,
                       RP = RP,
                       QTb = QTb,
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
} # End nnls_QR





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
  w <- crossprod(b, A) # calculate Gradient "t(A) %*% (b - A %*% x)" for vector x equal to zero 
  wmax <- max(w)
  Stellen <- -log(Fehlertol, 10)
  Mainloop <- 0

  #  InnerloopSum <- 0
  jprev <- 0
  Stabiler_Zustand <- 0
  Abbruch_Grund <- 0

  TestTol <- .Machine$double.eps * 1000
  if (Modellpunkte == 0) {
    Modellpunkte <- dim(A)[1]
  }
  ModellpunkteMax <- min(2*Modellpunkte, dim(A)[1])
  
  pb <- winProgressBar(title = "Fortschritt", min = 0, max = Modellpunkte, width = 600)
  
  # Mainloop
  while ((length(P) < ModellpunkteMax) & (Mainloop < MainloopMax) & (Stabiler_Zustand == 0) & (wmax > Fehlertol)) {
    Stabiltest <- 0
    while (Stabiltest == 0) {
      j <- which.max(w[R])
      Kandidat <- R[j] 
      Ptemp <- unique(sort(rbind(P, Kandidat)))
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
      sP <- solve( crossprod(AP), crossprod(AP, b))
#      sP <- solve( t(AP) %*% AP, t(AP) %*% b)
#      sP <- as.matrix(lsfit( AP, b, intercept = FALSE)[[1]])
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
                        Abbruch = Abbruch_Grund)
  
  class(nnls_VMF.out) <- "nnls"
  nnls_VMF.out 
} # End nnls_VMF





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





Givens_VMF <- function(Q, R, QTb, Spalte = 1) { # Givens Rotation
  Fehler <- 0
  
  if (Spalte > dim(R)[2]) {
    Fehler <- 1
  }
  
  if (Fehler == 0) {
    if (Spalte == (dim(R)[2])) { # Die letzte Spalte wird gestrichen => nichts zu berechnen
      RNeu <- R
      QNeu <- Q
      QTbNeu <- QTb
      print("Letze Spalte wird gestrichen")
    } else {
      n <- dim(R)[2]
      RNeu <- R
      QNeu <- Q
      QTbNeu <- QTb
      
      for (Lauf in (Spalte+1):n) { # Givens-Transformation
        r <- sign(RNeu[Lauf, Lauf]) * sqrt(RNeu[Lauf, Lauf]^2 + RNeu[Spalte,Lauf]^2)
        cc <- RNeu[Lauf, Lauf] / r
        ss <- RNeu[Spalte, Lauf] / r
        G <- matrix(c(cc,-ss,ss,cc),nrow = 2)
        
        v <- Lauf:n
        temp <- G %*% RNeu[c(Lauf,Spalte),v]
        RNeu[c(Lauf,Spalte),v] <- temp
        
        temp <- QNeu[,c(Lauf,Spalte)] %*% t(G)
        QNeu[,c(Lauf,Spalte)] <- temp
        
        temp <- G %*% QTbNeu[c(Lauf,Spalte)]
        QTbNeu[c(Lauf,Spalte)] <- temp
      } # Ende for-Schleife "Lauf
    } # Ende if "Spalte"-Abfrage
    
    Q <- QNeu[,-Spalte]
    R <- RNeu[-Spalte,-Spalte]
    QTb <- QTbNeu[-Spalte]
  } # Ende if "Fehler"-Abfrage    
  
  Givens_VMF.out <- list( Q = Q,
                          R = R,
                          QTb = QTb)
  Givens_VMF.out
}





Householder_VMF <- function(A, b, P, Spalte = 1) { # Householder Rotation
  Fehler <- 0
  
  if (Spalte > dim(A)[2]) {
    Fehler <- 1
  }
  
  if (sum(P == Spalte) > 0) {
    Fehler <- 1
  }
  
  if (Fehler == 0) {
    
    if (length(P) > 0) {
      R <- (1:N)[-P]
    } else {
      R <- (1:N)
    }
    
    P <- c(P, Spalte)
    n <- length(P)
    v <- dim(A)[1]
    N <- dim(A)[2]
    xneu <- A[n:v, Spalte]
    e <- as.matrix(c(1, rep(0, length(xneu)-1)))
    u <- sign(xneu[1]) * sqrt(sum(xneu^2)) * e + xneu
    
    A[n:v,R] <- A[n:v,R] - 2 * u %*% (t(u) %*% A[n:v,R]) / as.vector(t(u) %*% u)
    A[n:v, Spalte] <- round(A[n:v, Spalte],10)
    b[n:v] <- b[n:v] - 2 * u %*% (t(u) %*% b[n:v]) / as.vector(t(u) %*% u)
  } # Ende if "Fehler"-Abfrage    
  
  Householder_VMF <- list( A = A,
                           b = b,
                           P = P,
                           Fehler = Fehler)
  Householder_VMF
}





# Von der Seite https://www.r-bloggers.com/qr-decomposition-with-householder-reflections/
qr.householder <- function(A) {
  require(Matrix)
  
  R <- as.matrix(A) # Set R to the input matrix A
  
  n <- ncol(A)
  m <- nrow(A)
  H <- list() # Initialize a list to store the computed H matrices to calculate Q later
  V <- list() # Initialize a list to store the computed vk 
  QL <- list() # 
  RL <- list() # 
  
  if (m > n) {
    c <- n
  } else {
    c <- m
  }
  
  for (k in 1:c) {
    x <- R[k:m,k] # Equivalent to a_1
    e <- as.matrix(c(1, rep(0, length(x)-1)))
    vk <- sign(x[1]) * sqrt(sum(x^2)) * e + x
    
    # Compute the H matrix
    hk <- diag(length(x)) - 2 * vk %*% t(vk) / as.vector(t(vk) %*% vk)
    if (k > 1) {
      hk <- bdiag(diag(k-1), hk)
    }
    
    # Store the H matrix to find Q at the end of iteration
    H[[k]] <- hk
    V[[k]] <- vk
    
    R <- hk %*% R
    
    RL[[k]] <- R
  }
  
  Q <- Reduce("%*%", H) # Calculate Q matrix by multiplying all H matrices
  res <- list('Q'=Q,'R'=R,'RL'=RL,'H'=H,'V'=V)
  return(res)
} # End qr.householder


# Von der Seite "https://stackoverflow.com/questions/18349053/fastest-way-for-multiplying-a-matrix-to-a-vector!

func <- 'NumericMatrix mmult( NumericMatrix m , NumericVector v , bool byrow = true ){
  if( byrow );
    if( ! m.nrow() == v.size() ) stop("Non-conformable arrays") ;
  if( ! byrow );
    if( ! m.ncol() == v.size() ) stop("Non-conformable arrays") ;

  NumericMatrix out(m) ;

  if( byrow ){
    for (int j = 0; j < m.ncol(); j++) {
      for (int i = 0; i < m.nrow(); i++) {
        out(i,j) = m(i,j) * v[j];
      }
    }
  }
  if( ! byrow ){
    for (int i = 0; i < m.nrow(); i++) {
      for (int j = 0; j < m.ncol(); j++) {
        out(i,j) = m(i,j) * v[i];
      }
    }
  }
  return out ;
}'

