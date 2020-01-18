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

# source the scripts 
source("01_prepare_input.R")
source("02_prepare_output.R")
source("03_split_data_set.R")
source("04_apply_NN.R")


DT.input <- prepare_input()
DT.ouput <- prepare_output(DT.input$DGO_SPCODE_I)
DT.model <- split_data(DT.input, DT.ouput, 0.80)
learnNN(DT.model)

