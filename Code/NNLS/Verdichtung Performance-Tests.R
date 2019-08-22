#' ---
#' Title:   "The file contains some compression algorithm (Cluster, Optimisation)"
#' Author:  ""
#' Date:    "July 2019"
#' ---

start_time <- Sys.time()

#' ## General Settings

#' The pacman package is an R package management tool for loading and installing packages if necessary.
#'  
#' - **data.table:** For filtering, grouping and transforming the data as well as fast read and write 
#' - **dplyr:** For the piping operator.
#' - **lubridate:** Functions to work with date-times and time-spans
#' - **openxlsx:** High level interface to writing, reading and editing worksheets.
#' - **pivottabler:** Create regular pivot tables.
#' - **tcltk:** User-interface for selecting files and folders.
#' - **rstudioapi:** Setting the working directory to the source file location.
#' - **stringr:** Set of functions designed to make working with strings easy.
#' - **testthat:** Unit testing functions mainly used for sanity checks of the raw data.
#' - **crayon:** Color for terminal output, so that user defined errors appear in red. 
#' - **kableExtra:** Build common complex tables and manipulate table styles for rmarkdown. 
#' - **nnls:** The Lawson-Hanson algorithm for non-negative least squares (NNLS). 

#+ load-tes
suppressWarnings(library("pacman"))
pacman::p_load(data.table, dplyr, lubridate, openxlsx, pivottabler, 
               tcltk, rstudioapi, stringr, testthat, crayon, kableExtra, 
               nnls, matrixStats, microbenchmark, Rcpp, pracma, Rfast) #bit64

Hauptverzeichnis <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/")
setwd (Hauptverzeichnis)

source("Verdichtung Funktionen.R")



#
# Test Genauigkeit QR-Zerlegung
#

N <- 50000
v <- 500
n <- 100

A <- matrix(runif(N*v,0,1000000), nrow = v)
x <- rep(0,N)
P <- sort(sample(1:N,n))
R <- (1:N)[-P]
x[P] <- runif(n, 0, 1000)
b <- A %*% x + runif(v, 0, 0.5)

AP <- A[,P]
QP <- qr.Q(qr(AP))
RP <- qr.R(qr(AP))
max(abs(AP - QP %*% RP))

H <- householder(AP)
QH <- H$Q
RH <- H$R
max(abs(AP - QH %*% RH))

QTb <- t(QP) %*% b

x1 <- solve(crossprod(AP), crossprod(AP, b)) 
x2 <- matrix(backsolve(RP, t(QP) %*% b), ncol=1)
x3 <- matrix(backsolve(RH, t(QH) %*% b), ncol=1)
x4 <- matrix(backsolve(RP, QTb), ncol=1)
max(abs(x1-x2))
max(abs(x1-x3))
max(abs(x1-x4))
max(abs(x3-x4))

w1 <- t(A) %*% (b - AP %*% x1)
w2 <- t(A) %*% (b - AP %*% x2)
w3 <- t(A) %*% (b - AP %*% x3)
max(abs(w1-w2))
max(abs(w1-w3))
which.max(w1)
which.max(w2)

w1 <- t(A[,-P]) %*% (b - AP %*% x1)
w2 <- t(A[,-P]) %*% (b - AP %*% x2)
w3 <- t(A[,-P]) %*% (b - AP %*% x3)
max(abs(w1-w2))
max(abs(w1-w3))



#
# Performance-Tests "Least Square"-Berechnung
#

N <- 200
v <- 300
n <- 50

A <- matrix(runif(N*v,0,10000), nrow = v)
x <- rep(0,N)
P <- sort(sample(1:N,n))
AP <- A[,P]
x[P] <- runif(n, 0, 100)
b <- A %*% x

start.time <- Sys.time()
a1 <- solve( t(A) %*% A, t(A) %*% b)
end.time <- Sys.time()
end.time - start.time

start.time <- Sys.time()
AA <- t(A) %*% A
AB <- t(A) %*% b
a2 <- solve( AA, AB)
end.time <- Sys.time()
end.time - start.time

start.time <- Sys.time()
a3 <- as.matrix(lsfit( A, b, intercept = FALSE)[[1]])
end.time <- Sys.time()
end.time - start.time

start.time <- Sys.time()
A_QR <- qr(A, tol=1e-15)
Q <- qr.Q(A_QR)
R <- qr.R(A_QR)
a4 <- solve( R, t(Q) %*% b)
end.time <- Sys.time()
end.time - start.time

start.time <- Sys.time()
a5 <- solve(qr(A, LAPACK=TRUE), b)
end.time <- Sys.time()
end.time - start.time

start.time <- Sys.time()
a6 <- lm(b ~ 0 + A)
end.time <- Sys.time()
end.time - start.time

start.time <- Sys.time()
A_QR <- qr(A, tol=1e-15)
Q <- qr.Q(A_QR)
R <- qr.R(A_QR)
a7 <- backsolve( R, t(Q) %*% b)
end.time <- Sys.time()
end.time - start.time

microbenchmark(
  solve( t(A) %*% A ) %*% t(A) %*% b,
  solve( t(A) %*% A, t(A) %*% b),
  solve( crossprod(A), crossprod(A,b)),
  as.matrix(lsfit( A, b, intercept = FALSE)[[1]]),
  solve(qr(A, LAPACK=FALSE), b),
  solve(qr(A, LAPACK=TRUE), b),
  lm(b ~ 0 + A),
  solve( R, t(Q) %*% b),
  backsolve( R, t(Q) %*% b)
)

A1 <- t(A) %*% A
microbenchmark(
  solve( A1 ),
  spdinv( A1)
)

microbenchmark(
  t(A) %*% A,
  Rfast::mat.mult(Rfast::transpose(A), A)
)

microbenchmark(
  w <- t(A) %*% (b - A %*% x),
  w1 <- t(A) %*% (b - AP %*% sP),
  w2 <- crossprod(b - AP %*% sP, A)
)

A1 <- as(A, "dgeMatrix")
b1 <- as(b, "dgeMatrix")

microbenchmark(
  solve( t(A1) %*% A1 ) %*% t(A1) %*% b1,
  solve( t(A1) %*% A1, t(A1) %*% b1),
  solve( crossprod(A1), crossprod(A1,b1))
)

N <- 65000
n <- 300
x <- rep(0,N)
P <- sort(sample(1:N,n))



#
# Performance-Tests Vektor-Addition
#

N <- 50000
n <- 5000
x <- rep(0,N)
P <- sort(sample(1:N,n))
x[P] <- runif(n, 0, 100)

microbenchmark(
  x <- x + 0,
  x[P] <- x[P] + 0
)



#
# QR Remove Spalte - Performance-Tests
#

N <- 200
v <- 300
n <- 200

A <- matrix(runif(N*v,0,100), nrow = v)
PNeu <- 1:n

A_QR <- qr(A)
Q <- qr.Q(A_QR)
R <- qr.R(A_QR)

# Neuberechnung der gesamten QR-Zerlegung
f1 <- function(A, Q, R, Spalte = 1) { 
  ANeu <- A[,-Spalte]
  ANeu_QR <- qr(ANeu)
  QNeu <- qr.Q(ANeu_QR)
  RNeu <- qr.R(ANeu_QR)
  f1.out <- list( ANeu = ANeu,
                  QNeu = QNeu,
                  RNeu = RNeu)
}

# Rekursive Neuberechnung der QR-Zerlegung ab der Spalte die gestrichen wird
f2 <- function(A, Q, R, Spalte = 1) {
  ANeu <- A[,-Spalte]
  if (Spalte == 1) {
    RNeu <- NULL
    QNeu <- NULL
  } else {
    QNeu <- Q[,1:(Spalte-1)]
    RNeu <- R[1:(Spalte-1),1:(Spalte-1)]
  }
  
  for (Lauf in Spalte:(dim(ANeu)[2])) {
    if (Lauf == 1) {
      RNeu <- as.double(sqrt(t(ANeu[,1]) %*% ANeu[,1]))
      QNeu <- ANeu[,1] / RNeu
    } else {
      u <- ANeu[,Lauf]
      u <- u - QNeu %*% t(u %*% QNeu)
      unorm <- as.double(sqrt(t(u) %*% u))
      u <- u / unorm
      RNeu <- cbind(rbind(RNeu,0), c(ANeu[,Lauf] %*% QNeu, unorm))
      QNeu <- cbind(QNeu, u)
    }
  }
  
  f2.out <- list( ANeu = ANeu,
                  QNeu = QNeu,
                  RNeu = RNeu)
}

# Streichung mit Hilfe der "Givens Rotation"
f3 <- function(A, Q, R, Spalte = 1) { 
  n <- dim(A)[2]
  ANeu <- A[,-Spalte]
  RNeu <- R
  QNeu <- Q
  
  for (Lauf in (Spalte+1):n) {
    r <- sign(RNeu[Lauf, Lauf]) * sqrt(RNeu[Lauf, Lauf]^2 + RNeu[Spalte,Lauf]^2)
    cc <- RNeu[Lauf, Lauf] / r
    ss <- RNeu[Spalte, Lauf] / r
    G <- matrix(c(cc,-ss,ss,cc),nrow = 2)
    
    v <- Lauf:n
    temp <- G %*% RNeu[c(Lauf,Spalte),v]
    RNeu[c(Lauf,Spalte),v] <- temp
    
    temp <- QNeu[,c(Lauf,Spalte)] %*% t(G)
    QNeu[,c(Lauf,Spalte)] <- temp
  }
  
  QNeu <- QNeu[,-Spalte]
  RNeu <- RNeu[-Spalte,-Spalte]
  
  f3.out <- list( ANeu = ANeu,
                  QNeu = QNeu,
                  RNeu = RNeu)
}

Spalte <- 2
a <- f1 (A, Q, R, Spalte)
max(a$ANeu - a$QNeu %*% a$RNeu)

a <- f2 (A, Q, R, Spalte)
max(a$ANeu - a$QNeu %*% a$RNeu)

a <- f3 (A, Q, R, Spalte)
max(a$ANeu - a$QNeu %*% a$RNeu)

microbenchmark(
  f1 (A, Q, R, Spalte),
  f2 (A, Q, R, Spalte),
  f3 (A, Q, R, Spalte)
)

Spalte <- 150
a <- f3 (A, Q, R, Spalte)
max(a$ANeu - a$QNeu %*% a$RNeu)
Spalte <- 180
a <- f3 (a$A, a$Q, a$R, Spalte)
max(a$ANeu - a$QNeu %*% a$RNeu)



#
# Test - Berechnung der QR-Zerlegung
#

N <- 20
v <- 15
n <- 10
  
A <- matrix(runif(N*v,0,100), nrow = v)
x <- rep(0,N)
PNeu <- passive <- sort(sample(1:N,n))
i <- PNeu[n]
P <- PNeu[-n]
AP <- A[,P]
x[passive] <- runif(n, 0, 100)
b <- A %*% x
  
A_QR <- qr(A)
Q <- qr.Q(A_QR)
R <- qr.R(A_QR)

AP <- A[,P]
APNeu <- A[,PNeu]

AP_QR <- qr(AP)
QP <- qr.Q(AP_QR)
RP <- qr.R(AP_QR)

APNeu_QR <- qr(APNeu)
QPNeu <- qr.Q(APNeu_QR)
RPNeu <- qr.R(APNeu_QR)

u <- A[,i]
u <- u - QP %*% t(u %*% QP)
VZ <- -sign(u[length(PNeu)])
unorm <- VZ * as.double(sqrt(t(u) %*% u))
u <- u / unorm

RPNeuber <- cbind(rbind(RP,0), c(A[,i] %*% QP, unorm))
QPNeuber <- cbind(QP, u)

max(abs(QPNeu - QPNeuber))
max(abs(RPNeu - RPNeuber))

u[v]
u[n]

QPNeu[v,n]
QPNeuber[v,n]

RPNeu[n,n]
RPNeuber[n,n]

Stellen <- 10
max(round(AP - QP %*% RP, Stellen))
max(round(APNeu - QPNeu %*% RPNeu, Stellen))
max(round(APNeu - QPNeuber %*% RPNeuber, Stellen))

# QR update
for (Lauf in 1:length(PNeu)) {
  if (Lauf == 1) {
    APLauf <- A[,PNeu[Lauf]]
    RPLauf <- as.double(sqrt(t(APLauf) %*% APLauf))
    QPLauf <- APLauf / RPLauf
    QTb <- QPLauf %*% b
  } else {
    u <- A[,PNeu[Lauf]]
    u <- u - QPLauf %*% t(u %*% QPLauf)
    unorm <- as.double(sqrt(t(u) %*% u))
    u <- u / unorm
    RPLauf <- cbind(rbind(RPLauf,0), c(A[,PNeu[Lauf]] %*% QPLauf, unorm))
    QPLauf <- cbind(QPLauf, u)
    QTb <- c(QTb, t(u) %*% b)
  }
}

Stellen <- 4
round(QPNeu, Stellen)
round(QPLauf, Stellen)

round(RPNeu, Stellen)
round(RPLauf, Stellen)

round(QPNeu %*% RPNeu, Stellen) - round(QPLauf %*% RPLauf, Stellen)

t(QPLauf) %*% b
QTb

sP <- backsolve(RPLauf, t(QPLauf) %*% b)
s <- rep(0, N)
s[PNeu] <- sP

# QR Remove Spalte => Test Givens-Transformation
Spalte <- 5
n <- length(PNeu)
R <- RPNeu
Q <- QPNeu
for (Lauf in (Spalte+1):n) {
  r <- sign(R[Lauf, Lauf]) * sqrt(R[Lauf, Lauf]^2 + R[Spalte,Lauf]^2)
  cc <- R[Lauf, Lauf] / r
  ss <- R[Spalte, Lauf] / r
  G <- matrix(c(cc,-ss,ss,cc),nrow = 2)

  v <- Lauf:n
  temp <- round(G %*% R[c(Lauf,Spalte),v],8)
  R[c(Lauf,Spalte),v] <- temp
  
  temp <- round(Q[,c(Lauf,Spalte)] %*% t(G),8)
  Q[,c(Lauf,Spalte)] <- temp
}

Q[,-Spalte] %*% R[-Spalte, -Spalte]
APNeu[,-Spalte]
max(abs(Q[,-Spalte] %*% R[-Spalte, -Spalte] - APNeu[,-Spalte]))

R[-Spalte,-Spalte]
qr.R(qr(A[,PNeu[-Spalte]]))
max(abs(R[-Spalte,-Spalte]) - abs(qr.R(qr(A[,PNeu[-Spalte]]))))

Stellen <- 4
round(qr.Q(qr(A[,PNeu])), Stellen)
round(qr.Q(qr(A[,PNeu[-Spalte]])), Stellen)

round(qr.R(qr(A[,PNeu])), Stellen)
round(qr.R(qr(A[,PNeu[-Spalte]])), Stellen)



#
# Test Matrixenmultiplikation
#

N <- 500
n <- 50

A <- matrix(runif(N*n,0,100), nrow = n)
b <- runif(n,0,100)

microbenchmark(
  crossprod(A),
  t(A) %*% A
)

microbenchmark(
  crossprod(A, b),
  t(A) %*% b
)



#
# Test Gradientenberechnung
#

N <- 10000
v <- 500
n <- 100

A <- matrix(runif(N*v,0,100), nrow = v)
x <- rep(0,N)
P <- sort(sample(1:N,n))
AP <- A[,P]
x[P] <- runif(n, 0, 100)
b <- A %*% x + runif(v, 0, 10)
sP <- x[P]

Gradfix <- crossprod(b, A)

# Berechnung der Gradienten gesamt
f1 <- function(A, AP, b, sP) {
  crossprod(b - AP %*% sP, A)
}

# Berechnung des Gradienten Stückweise
f2 <- function(A, AP, b, sP) {
  temp <- AP %*% sP
  temp2 <- crossprod(temp, A)
  Gradfix - temp2
}

microbenchmark(
  f1(A, AP, b, sP),
  f2(A, AP, b, sP)
)

microbenchmark(
  crossprod(b - AP %*% sP, A),
  t(b - AP %*% sP) %*% A,
  Gradfix - t(AP %*% sP) %*% A
)






# Test Householder

source("Verdichtung Funktionen.R")

library (pracma)

N <- 20
v <- 15
n <- 5

A <- matrix(runif(N*v,0,100), nrow = v)
x <- rep(0,N)
Porig <- sort(sample(1:N,n))
AP <- A[,Porig]
x[Porig] <- runif(n, 0, 100)
b <- A %*% x + runif(v, 0, 10)

Erg <- nnls_Fortran(A, b)
Porig <- Erg$passive

ASicher <- A
bSicher <- b

A <- ASicher
b <- bSicher
v <- dim(A)[1]
P <- NULL

for (Lauf in 1:length(Porig)) {
  i <- Porig[Lauf]
  P <- c(P, i)
  R <- (1:N)[-P]
  n <- length(P)
  x <- A[n:v,i]
  e <- as.matrix(c(1, rep(0, length(x)-1)))
  u <- sign(x[1]) * sqrt(sum(x^2)) * e + x
  
  # Compute the H matrix
  H <- eye(length(x)) - 2 * u %*% t(u) / as.vector(t(u) %*% u)
  if (n > 1) {
    H <- rbind(cbind(eye(n-1),zeros(n-1,length(x))),cbind(zeros(length(x),n-1),H))
  }
  A <- H %*% A
  b <- t(H) %*% b
}  

#R <- A[1:length(P),P]
#round(matrix(Erg$A, nrow=v)[1:Erg$nsetp,Erg$passive])

Aalt <- A
balt <- b
Palt <- P 
i <- 2 # Spalte die gelöscht werden soll

P <- Palt
A <- Aalt
b <- balt
Ralt <- round(A[1:length(P),P],5)
Ralt

for (Lauf in i:(n-1)) {
  r <- sign(A[Lauf, P[Lauf]]) * sqrt(A[Lauf, P[Lauf]]^2 + A[Lauf+1,P[Lauf]]^2)
  cc <- A[Lauf, P[Lauf]] / r
  ss <- A[Lauf+1, P[Lauf]] / r
  G <- matrix(c(cc,-ss,ss,cc),nrow = 2)
  
  temp <- G %*% A[c(Lauf,Lauf+1),]
  A[c(Lauf,Lauf+1),] <- temp

  temp <- G %*% b[c(Lauf,Lauf+1)]
  b[c(Lauf,Lauf+1)] <- temp
  
}
round(A[1:length(P),P],5)

P <- P[-i]



Ralt
round(A[1:length(P),P],5)

P <- P[-i]
Ralt
round(A[1:length(P)+1,P],5)


A <- ASicher
b <- bSicher
v <- dim(A)[1]
P <- NULL

for (Lauf in 1:length(Porig)) {
  i <- Porig[Lauf]
  a <- Householder_VMF(A,b,P,i)
  P <- c(P, i)
  A <- a$A
  b <- a$b
  P <- a$P
}  

A[,P]



#
# Performance Householder-Transformation
# Multiplikation H * A
#

N <- 200
v <- 50
i <- 3
A <- matrix(runif(N*v,0,100), nrow = v)
b <- matrix(runif(v,0,100), nrow = v)
n <- 1
x <- A[n:v,i]
e <- as.matrix(c(1, rep(0, length(x)-1)))
u <- sign(x[1]) * sqrt(sum(x^2)) * e + x

# Compute the H matrix
H <- eye(length(x)) - 2 * u %*% t(u) / as.vector(t(u) %*% u)
if (n > 1) {
  H <- rbind(cbind(eye(n-1),zeros(n-1,length(x))),cbind(zeros(length(x),n-1),H))
}

microbenchmark(
  H %*% A,
  A - 2 * u %*% (t(u) %*% A) / as.vector(t(u) %*% u)
)

microbenchmark(
  t(H) %*% b,
  b - 2 * u %*% (t(u) %*% b) / as.vector(t(u) %*% u)
)

t(H) %*% b
b - 2 * u %*% (t(u) %*% b) / as.vector(t(u) %*% u)



PNeu <- passive <- sort(sample(1:N,n))
P <- PNeu[-n]
i <- PNeu[n]
APNeu <- A[,PNeu]
AP <- A[,P]

APNeu_QR <- qr(APNeu)
QPNeu <- qr.Q(APNeu_QR)
RPNeu <- qr.R(APNeu_QR)

AP_QR <- qr(AP)
QP <- qr.Q(AP_QR)
RP <- qr.R(AP_QR)

x <- A[n:v,i]
e <- as.matrix(c(1, rep(0, length(x)-1)))
u <- sign(x[1]) * sqrt(sum(x^2)) * e + x

# Compute the H matrix
H <- eye(length(x)) - 2 * u %*% t(u) / as.vector(t(u) %*% u)
H <- rbind(cbind(eye(n-1),zeros(n-1,length(x))),cbind(zeros(length(x),n-1),H))

cbind(QP,A[,i])

a <- qr.householder(APNeu)

