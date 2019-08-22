#' ---
#' Title:   "The file contains some compression algorithm (Cluster, Optimisation)"
#' Author:  ""
#' Date:    "July 2019"
#' ---

start_time <- Sys.time()

#+ global_options, include = FALSE
# knitr::opts_chunk$set(eval = FALSE)

# Set digits
# options("digits" = 4)

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

#cppFunction( func )

library(parallel)
cores <- detectCores() 
print(cores)

# 
# Analysis for "IF_G_DGO"
#

# Read in variablen

Steuerung_Cluster_Variablen <- read.xlsx(xlsxFile = paste0(Hauptverzeichnis, "Steuerung Cluster.xlsx"), 
                               namedRegion = "Variablen")
Steuerung_Cluster_Auswertungsjahre <- read.xlsx(xlsxFile = paste0(Hauptverzeichnis, "Steuerung Cluster.xlsx"),
                                                namedRegion = "Auswertungsjahre")
Variablen <- Steuerung_Cluster_Variablen %>% 
  filter (Verwenden == 1) %>% .$Variable # Contains all variablen which should be included 
Standardisieren <- Steuerung_Cluster_Variablen %>% 
  filter (Verwenden == 1, Standardisieren == 1) %>% .$Variable
Auswertungsjahre <- as.integer(Steuerung_Cluster_Auswertungsjahre$Jahr)
rm (Steuerung_Cluster_Variablen, Steuerung_Cluster_Auswertungsjahre)

#
# Read in data
#

# Policy-informations
OrigDaten <- paste0(Hauptverzeichnis, "IF_G_DGO/") %>% 
  paste0(., "DMESGP.txt") %>%
  fread(file = ., dec = ".", sep = ",", skip = 0,  stringsAsFactors = TRUE, blank.lines.skip = TRUE) %>%
  setDT()

DatenPro <- paste0(Hauptverzeichnis, "RUN_986/DMESGP/") %>%
  paste0(., "DMESGP~ACT_W_IF_DMESGP~Base.gro") %>%
  fread(file = ., dec = ".", sep = ",", skip = 2,  stringsAsFactors = TRUE, blank.lines.skip = TRUE) %>%
  select(Variablen) %>% 
  setDT()

Sollwerte <- paste0(Hauptverzeichnis, "RUN_984/DMESGP/") %>%
  paste0(., "DMESGP~ACT_W_IF_DMESGP~Base.ref") %>%
  fread(file = ., dec = ".", sep = ",", skip = 2,  stringsAsFactors = TRUE, blank.lines.skip = TRUE) %>%
  select(Variablen) %>% 
  setDT()

# Numbers of SPCodes
MaxSPCode <- max(DatenPro[,DGO_SPCODE_I])

# Add the year 
#Jahr <- unlist(lapply(diff(c(NeuSPCode,dim(DatenPro)[1]+1)),seq)) - 1 # Old version
DatenPro[, Jahr := 1:.N - 1, by = DGO_SPCODE_I]
Sollwerte[, Jahr := 1:.N - 1, by = DGO_SPCODE_I]

# Determine the lines in which a SP code appears for the first time
NeuSPCode <- match(1:MaxSPCode, DatenPro[,DGO_SPCODE_I]) # no longer needed



#
# Some tests => can be removed
#

OrigDaten[,.(AnzPol = sum(INIT_POLS_IF_I), 
             AnzStamm = sum(STAMM_SEG_I * INIT_POLS_IF_I),
             BilRes = sum(BIL_RES_I * INIT_POLS_IF_I), 
             Pämie = sum(ANN_PRMAL_PP_I * INIT_POLS_IF_I),
             EE = sum(SINGPRMAL_PP_I * INIT_POLS_IF_I),
             VAGRes = sum(FONDSW_VAG_I * INIT_POLS_IF_I),
             Fonds = sum(INIT_FNDWERT_I * INIT_POLS_IF_I),
             Garantie = sum(GARANTIE_I * INIT_POLS_IF_I))]

DatenPro [Jahr == 0, .(AnzPol = sum(NO_POLS_IF), 
                       AnzStamm = sum(NO_STAMM_IF),
                       BilRes = sum(MATH_RES_IF), 
                       Pämie = sum(PREM_INC),
                       VAGRes = sum(C_MATHRES_IF),
                       Fonds = sum(U_MATHRES_IF))]

Sollwerte [Jahr == 0, .(AnzPol = sum(NO_POLS_IF), 
                       AnzStamm = sum(NO_STAMM_IF),
                       BilRes = sum(MATH_RES_IF), 
                       Pämie = sum(PREM_INC),
                       VAGRes = sum(C_MATHRES_IF),
                       Fonds = sum(U_MATHRES_IF))]


#StandDatenPro <- copy(DatenPro)
#StandDatenPro[, (Standardisieren) := lapply(.SD, scale), .SDcols=Standardisieren]

Spalten <- colnames(DatenPro)[c(-1, -dim(DatenPro)[2])]
Daten <- dcast(DatenPro[Jahr %in% Auswertungsjahre, ], DGO_SPCODE_I ~ Jahr, value.var = Spalten) %>% 
  select(., -DGO_SPCODE_I) %>% as.matrix(.) %>% t(.)
Daten[is.na(Daten)] <- 0 # SPCodes mit geringerer Laufzeit als das Maximum in "Auswertungsjahre" haben "na" als Wert.

SollDaten <- dcast(Sollwerte[Jahr %in% Auswertungsjahre, ], DGO_SPCODE_I ~ Jahr, value.var = Spalten) %>%
  select(., -DGO_SPCODE_I) %>% as.numeric(.) %>% t(.)


# Normalize "Daten" and "SollDaten"

if (1 == 2) {
RMax <- rowMaxs(Daten)
RMin <- rowMins(Daten)
Normierung <- RMax != RMin

Daten[Normierung, ] <- (Daten[Normierung, ] - RMin[Normierung]) / (RMax[Normierung] - RMin[Normierung])
SollDaten[Normierung] <- (SollDaten[Normierung] - RMin[Normierung]) / (RMax[Normierung] - RMin[Normierung])

rm(RMax, RMin, Normierung)
}

DatenSicherung <- Daten

Daten <- DatenSicherung
source("Verdichtung Funktionen.R")

N <- 20000
#N <- dim(Daten)[2]
Daten <- Daten[,1:N]

system.time(Erg <- nnls(Daten, t(SollDaten)))
Erg$nsetp
sort(Erg$passive)

# Input

A <- Daten
b <- t(SollDaten)

system.time(a <- nnls_QR(A, b, 0, Ausgabe = 0, MainloopMax = 0, Test = 1))
a$nsetp
sort(a$passive)

system.time(a1 <- nnls_VMF(A, b, 0, Ausgabe = 0))
a1$nsetp
sort(a1$passive)



system.time(a2 <- nnls_QR_House(A, b, 0, Ausgabe = 0, MainloopMax = 0, Test = 1))
a2$nsetp
sort(a2$passive)

