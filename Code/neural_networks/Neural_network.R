############## load all required packages ##############

#Loading libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rstudioapi, data.table, tidyverse, caret, ggplot2,
               viridis, gridExtra, pROC, descr, keras, mltools)


############## Set up the workspace for reproduceability ##############

# set the seed to get similar results for the random-parts
set.seed(100) 

# set the working directory to the path of the script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

############## Read and prepare the input-data ##############
# read the modelpoints
raw_input <- fread("EB/DMWV__.txt", header = TRUE, skip = "!")

# preselct those columns which could be important 
input_selCol <- raw_input[, .(DGO_SPCODE_I,
                              AGE_AT_ENTRY_I, B2_PRAEMIE_I, BIL_RES_I,
                              DURATIONIF_M_I, INIT_ACC_BON_I, INIT_POLSTAT_I,
                              POL_TERM_Y_I, PREM_FREQ_I, PREM_PAYBL_Y_I,
                              RABATT_KZ_I, SEX_I, STAMM_SEG_I, SUM_ASSURED_I,
                              VERTRIEB_I, ZUGANG_ART_I, ZV_PRAEMIE1_I)]

# one hot encoding for SEX_I
oneHot_SEX_I <- input_selCol[, .(DGO_SPCODE_I, SEX_I)]
oneHot_SEX_I[, SEX_I := as.character(SEX_I)]
oneHot_SEX_I[SEX_I == "0", SEX_I := "male"]
oneHot_SEX_I[SEX_I == "1", SEX_I := "female"]
oneHot_SEX_I_final <- dcast(data = oneHot_SEX_I,
                            DGO_SPCODE_I ~ SEX_I,
                            fun = length)

# one hot encoding for STAMM_SEG_I
oneHot_STAMM_SEG_I <- input_selCol[, .(DGO_SPCODE_I, STAMM_SEG_I)]
oneHot_STAMM_SEG_I[, STAMM_SEG_I := as.character(STAMM_SEG_I)]
oneHot_STAMM_SEG_I[STAMM_SEG_I == "0", STAMM_SEG_I := "NoStammSeg"]
oneHot_STAMM_SEG_I[STAMM_SEG_I == "1", STAMM_SEG_I := "StammSeg"]
oneHot_STAMM_SEG_I_final <- dcast(data = oneHot_STAMM_SEG_I,
                                  DGO_SPCODE_I ~ STAMM_SEG_I,
                                  fun = length)



# one hot encoding for INIT_POLSTAT
oneHot_INIT_POLSTAT_I <- input_selCol[, .(DGO_SPCODE_I, INIT_POLSTAT_I)]
oneHot_INIT_POLSTAT_I <- oneHot_INIT_POLSTAT_I[, INIT_POLSTAT_I := as.character(INIT_POLSTAT_I)]
oneHot_INIT_POLSTAT_I[INIT_POLSTAT_I == "1",  INIT_POLSTAT_I := "Polpämpfl"]
oneHot_INIT_POLSTAT_I[INIT_POLSTAT_I == "3",  INIT_POLSTAT_I := "Polpämfr"]
oneHot_INIT_POLSTAT_I[INIT_POLSTAT_I == "4",  INIT_POLSTAT_I := "Polpämfr2"]

oneHot_INIT_POLSTAT_I_final <- dcast(data = oneHot_INIT_POLSTAT_I,
                                     DGO_SPCODE_I ~ INIT_POLSTAT_I,
                                     fun = length)

# one hot encoding for PREM_FREQ_I --> nur Jahreszahler? 
oneHot_PREM_FREQ_I <- input_selCol[, .(DGO_SPCODE_I, PREM_FREQ_I)]
oneHot_PREM_FREQ_I <- oneHot_PREM_FREQ_I[, PREM_FREQ_I := as.character(PREM_FREQ_I)]
oneHot_PREM_FREQ_I[PREM_FREQ_I == "1",  PREM_FREQ_I := "jährlich"]
oneHot_PREM_FREQ_I[PREM_FREQ_I == "12", PREM_FREQ_I := "mtl"]
oneHot_PREM_FREQ_I[PREM_FREQ_I == "2",  PREM_FREQ_I := "halbjäh"]
oneHot_PREM_FREQ_I[PREM_FREQ_I == "4",  PREM_FREQ_I := "quartä"]

oneHot_PREM_FREQ_I_final <- dcast(data = oneHot_PREM_FREQ_I,
                                  DGO_SPCODE_I ~ PREM_FREQ_I,
                                  fun = length)

# one hot encoding for RABATT_KZ_I --> nur Fälle mit "N" nehmen! --> kein Encoding
input[, .N, by = .(RABATT_KZ_I)]

# one hot encoding for VERTRIEB_I --> Nur MAKLER
oneHot_VERTRIEB_I <- input_selCol[, .(DGO_SPCODE_I, VERTRIEB_I)]

oneHot_VERTRIEB_I_final <- dcast(data = oneHot_VERTRIEB_I,
                                  DGO_SPCODE_I ~ VERTRIEB_I,
                                  fun = length)

# one hot encoding for ZUGANG_ART_I --> in Prophet alles gleich behandelt solange
#                                       K oder P
input_selCol[, .N, by = .(ZUGANG_ART_I)]

#" 0" = Neuzugang
#" 1" = Anpassungsbrief (S) / Dynamisierung (D)
#" 2" = Ersatz (mit Rückdatierung)
#" 4" = Reduktion
#" 5" = Reaktivierung
#" 8" = Indexklauselerhöhung bei Großleben Anpassung bei Kollektiv (S)
#" 9" = Interner Zugang (infolge Ablauf der Prämienzahlung)
#" K" = Konvertierung (S) / Ersatz bei prämienfreier Summe (D)
#" N" = Nachtrag, Aufstockung, Indexbrief
#" P" = Prämienfreie Einrechnung (S)
#" S" = Prolongation (nur Großleben) (S) / Verlängerung (D)
#" T" = Versorgertod
#" W" = Switch ohne Tarifänderung ( FLV )
#" Y" = Switch mit Tarifänderung ( FLV )
#" Z" = Summenabfall wegen Zwischenauszahlung
#"  " bei Änderungsarten, die keinen Zugang inkludieren


#input_selCol[, .N,
#      by = .(INIT_POLSTAT_I, PREM_FREQ_I, RABATT_KZ_I, VERTRIEB_I, ZUGANG_ART_I)] %>% fwrite("analyse.csv") 

input_tmp <- cbind(input_selCol,
                   oneHot_SEX_I_final[, .(female, male)],
                   oneHot_STAMM_SEG_I_final[, .(NoStammSeg, StammSeg)])


# take a sample which fits our needs and doesnt need one hot encoding except for sex and STAMM_SEG_I
input_tmp <- input_tmp[INIT_POLSTAT_I == 1
                       & PREM_FREQ_I    == 12
                       & RABATT_KZ_I    == "N"
                       & VERTRIEB_I     == "MAKLER"
                       & ZUGANG_ART_I %in% c("N", "0"),]

input <- input_tmp[, .(DGO_SPCODE_I,
                       AGE_AT_ENTRY_I,
                       B2_PRAEMIE_I,
                       BIL_RES_I,
                       DURATIONIF_M_I,
                       INIT_ACC_BON_I,
                       POL_TERM_Y_I,
                       PREM_PAYBL_Y_I,
                       male,        # one hot from SEX_I
                       female,      # one hot from SEX_I
                       StammSeg,    # one hot from STAMM_SEG_I
                       NoStammSeg,  # one hot from STAMM_SEG_I
                       SUM_ASSURED_I,
                       ZV_PRAEMIE1_I)]



lapply(input, summary)






############## Read and prepare the output-data ##############
# read the raw-data output data
DT.output <- fread("RUN_986/D__RIS~ACT_W_EB_D__RIS~Base.gro",
                   sep = ",",
                   header = TRUE,
                   skip = "DGO_SPCODE_I")

# select and MATH_RES_IF as output and DGO_SPCODE_I as policy identifier
DT.output_sel <- DT.output[, .(DGO_SPCODE_I, MATH_RES_IF)]
# select only those rows corresponding to the input
DT.output_sel <- DT.output_sel[DGO_SPCODE_I %in% input$DGO_SPCODE_I, ]
# add increasing Jahr-column per DGO_SPCODE 
DT.output_sel[, Jahr := 0:.N, by = .(DGO_SPCODE_I)]
# output-matrix
DT.output_mat <- dcast(DT.output_sel, DGO_SPCODE_I ~ Jahr , value.var = "MATH_RES_IF", fill = 0)



#####################  Split data #######################

# Find randomly selected rownumbers of the raw data so that we get 90% of data for training
train_rows <- createDataPartition(DT.rawData$XXX, p = 0.9, list = FALSE)
# create list of data.tables for training purposes
DT.train <- apply(train_rows, 2, function(x){DT.rawData[x,]})
# create list of data.tables for testing purposes
DT.test <- apply(train_rows, 2, function(x){DT.rawData[-x,]})


learnNN <- function(DT.train, DT.test, columns, ...){
  
  x_train <- DT.train$Resample1[, ..columns]
  x_train <- x_train[ , lapply(.SD, log1p), .SDcols = 1:ncol(x_train)]
  x_train <- x_train[ , lapply(.SD, scale), .SDcols = 1:ncol(x_train)]
  
  # initialize the NN-model
  model_keras <- keras_model_sequential()
  
  # set up and compile the NN-model with the following parameters:
  model_keras %>% 
    # First hidden layer
    layer_dense(
      units              = 16, 
      kernel_initializer = "uniform", 
      activation         = "relu", 
      input_shape        = ncol(x_train)) %>% 
    
    # Dropout to prevent overfitting
    layer_dropout(rate = 0.1) %>%
    
    # Second hidden layer
    layer_dense(
      units              = 16, 
      kernel_initializer = "uniform", 
      activation         = "relu") %>% 
    
    # Dropout to prevent overfitting
    layer_dropout(rate = 0.1) %>%
    
    # Output layer
    layer_dense(
      units              = 1, 
      kernel_initializer = "uniform", 
      activation         = "sigmoid") %>% 
    
    # Compile ANN
    compile(
      optimizer = 'adam',
      loss      = 'binary_crossentropy',
      metrics   = c('accuracy')
    )
  
  # get 0 and 1 for non-spam and spam
  y_train <- as.integer(DT.train$Resample1[, is_spam_flag]) - 1
  
  trained_model <- fit(
    object = model_keras,
    x = as.matrix(x_train), 
    y = y_train, 
    epochs = 35, 
    batch_size = 32
  )
  model_keras
  x_test <- DT.test$Resample1[, ..columns]
  x_test <- x_test[ , lapply(.SD, log1p), .SDcols = 1:ncol(x_test)]
  x_test <- x_test[ , lapply(.SD, scale), .SDcols = 1:ncol(x_test)]
  
  
  
  erg_prob <- predict_proba(object = model_keras, x = as.matrix(x_test))
  erg_class <- predict_classes(object = model_keras, x = as.matrix(x_test))
  roc_model <- roc(response = DT.test$Resample1$is_spam_flag,
                   predictor = as.numeric(erg_prob))
  
  error <- sum(erg_class == (as.numeric(DT.test$Resample1$XXX) - 1)) / length(erg_class)
  
  return(list("Error Rate" = error, "roc" = roc_model))
}