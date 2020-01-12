############## load all required packages ##############

#Loading libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rstudioapi, data.table, tidyverse, caret, ggplot2,
               viridis, gridExtra, pROC, descr, )


############## Read and prepare the raw-data ##############

# set the seed to get similar results for the random-parts
set.seed(100) 

# set the working directory to the path of the script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read the header of the data
DT.header <- fread("XXX/XXX.txt", header = FALSE)
# read the raw-data
DT.rawData <- fread("XXX/XXX.txt", sep = ",", header = FALSE)



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