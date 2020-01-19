############## load all required packages ##############

#Loading libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rstudioapi, data.table, tidyverse, caret, ggplot2,
               viridis, gridExtra, pROC, descr, keras, mltools)


############## Set up the workspace for reproduceability ##############

# set the seed to get similar results for the random-parts
set.seed(120) 

# set the working directory to the path of the script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# source the scripts 
source("01_prepare_input.R")
source("02_prepare_output.R")
source("03_split_data_set.R")
source("04_apply_NN_01.R")
source("04_apply_NN_02.R")
source("04_apply_NN_03.R")
source("04_apply_NN_04.R")
source("04_apply_NN_05.R")
source("04_apply_NN_06.R")
source("04_apply_NN_07.R")
source("04_apply_NN_08.R")
source("05_plot_data.R")


DT.input <- prepare_input()
DT.ouput <- prepare_output(DT.input$DGO_SPCODE_I)
DT.model <- split_data(DT.input, DT.ouput, 0.80)

############### Three hidden layers ############### 

#flag_integer("dense_units1", 128)

result_01 <- lapply(1:40, function(x) learnNN_01(DT.model, 8, 'adam'))
#plot_data(result_01, "result_01.pdf")

result_02 <- lapply(1:40, function(x) learnNN_01(DT.model, 20, 'adam'))
plot_data(result_02, "result_02.pdf")

result_03 <- lapply(1:40, function(x) learnNN_01(DT.model, 37, 'adam'))
plot_data(result_03, "result_03.pdf")

result_04 <- lapply(1:40, function(x) learnNN_01(DT.model, 50, 'adam'))
plot_data(result_04, "result_04.pdf")

############### Two hidden layers ###############

result_05 <- lapply(1:40, function(x) learnNN_02(DT.model, 8, 'adam'))
plot_data(result_05, "result_05.pdf")

result_06 <- lapply(1:40, function(x) learnNN_02(DT.model, 20, 'adam'))
plot_data(result_06, "result_06.pdf")

result_07 <- lapply(1:40, function(x) learnNN_02(DT.model, 37, 'adam'))
plot_data(result_07, "result_07.pdf")

result_08 <- lapply(1:40, function(x) learnNN_02(DT.model, 50, 'adam'))
plot_data(result_08, "result_08.pdf")


