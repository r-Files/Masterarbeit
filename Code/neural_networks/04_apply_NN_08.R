
############## Set up neural network and test the results  ##############

learnNN_08 <- function(data, ...){
  
  # extract the train information
  x_train <- data$train$input
  y_train <- data$train$output
  
  # extract the test information  
  x_test <- data$test$input
  y_test <- data$test$output

  set.seed(120)
  
  cols <- c("male","female","StammSeg","NoStammSeg")

  # scale everything from input except the one hot encoded columns and then add those columns
  
  # training set
  x_train <- x_train[ , lapply(.SD, scale), .SDcols = -cols] %>% 
    cbind(data$train$input[, ..cols])
  
  # test set
  x_test <- x_test[ , lapply(.SD, scale), .SDcols = -cols] %>% 
    cbind(data$test$input[, ..cols])
  
  
  # initialize the NN-model
  model_keras <- keras_model_sequential()
  
  # set up and compile the NN-model with the following parameters:
  model_keras %>% 
    # First hidden layer
    layer_dense(
      units              = 50, 
      kernel_initializer = "uniform", 
      activation         = "relu", 
      input_shape        = ncol(x_train)) %>% 
    
    # Dropout to prevent overfitting
    layer_dropout(rate = 0.1) %>%
    
    # Second hidden layer
    layer_dense(
      units              = 50, 
      kernel_initializer = "uniform", 
      activation         = "relu") %>% 
    
    # Dropout to prevent overfitting
    layer_dropout(rate = 0.1) %>%
    
    # Output layer
    layer_dense(
      units              = ncol(y_train), 
      kernel_initializer = "uniform", 
      activation         = "relu") %>% 
    
    # Compile ANN
    compile(
      optimizer = 'adam',
      loss      = 'mean_squared_error',
      metrics   = c('mse')
    )
  
  trained_model <- fit(
    object = model_keras,
    x = as.matrix(x_train), 
    y = as.matrix(y_train), 
    epochs = 50, 
    batch_size = 32
  )
  
  # print info of model
  # model_keras

  result <- predict(object = model_keras, as.matrix(x_test))
  
  data.table(
    Year = 0:(length(colSums(result)) - 1),
    predicted = colSums(result),
    true = colSums(y_test)
  ) %>% return()
  
}