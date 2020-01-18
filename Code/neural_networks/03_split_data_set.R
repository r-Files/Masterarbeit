
############## Split input and output data for training and testing  ##############

split_data <- function(input, output, perc) {
  
  # determine the number of policies in the data-set 
  nDataSet <- nrow(input)
  
  # Find randomly selected rownumbers: use (perc)% of the data for training and 
  # (1 - perc)% for testing
  train_rows <- sample(1:nDataSet, nDataSet * perc)
  
  # create data.tables for training purposes and drop the column of the unique policy ID's
  DT.train_input  <- input[train_rows, !"DGO_SPCODE_I"]
  DT.train_output <- output[train_rows, !"DGO_SPCODE_I"]
  train <- list(input = DT.train_input,
                output = DT.train_output)
  
  # create data.tables for testing purposes
  DT.test_input  <- input[-train_rows, !"DGO_SPCODE_I"]
  DT.test_output <- output[-train_rows, !"DGO_SPCODE_I"]
  test <- list(input = DT.test_input,
               output = DT.test_output)
  
  return(list(train = train, test = test))

}
