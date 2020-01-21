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
source("06_analyse_results.R")


DT.input <- prepare_input()
DT.ouput <- prepare_output(DT.input$DGO_SPCODE_I)
DT.model <- split_data(DT.input, DT.ouput, 0.80)

############### Three hidden layers ############### 

set.seed(100)
result_01 <- lapply(1:40, function(x) learnNN_01(DT.model, 8, 'adam'))
saveRDS(result_01, file = "result_01.rds")
analyse_01 <- analyse_nn(result_01)
#plot_data(result_01, "../../figures/chapter_NN/result_01.pdf")

set.seed(100)
result_02 <- lapply(1:40, function(x) learnNN_01(DT.model, 20, 'adam'))
saveRDS(result_02, file = "result_02.rds")
analyse_02 <- analyse_nn(result_02)
#plot_data(result_02, "../../figures/chapter_NN/result_02.pdf")

set.seed(100)
result_03 <- lapply(1:40, function(x) learnNN_01(DT.model, 37, 'adam'))
saveRDS(result_03, file = "result_03.rds")
analyse_03 <- analyse_nn(result_03)
#plot_data(result_03, "../../figures/chapter_NN/result_03.pdf")

set.seed(100)
result_04 <- lapply(1:40, function(x) learnNN_01(DT.model, 50, 'adam'))
saveRDS(result_04, file = "result_04.rds")
analyse_04 <- analyse_nn(result_04)
#plot_data(result_04, "../../figures/chapter_NN/result_04.pdf")

############### Two hidden layers ###############
set.seed(100)
result_05 <- lapply(1:40, function(x) learnNN_02(DT.model, 8, 'adam'))
saveRDS(result_05, file = "result_05.rds")
analyse_05 <- analyse_nn(result_05)
#plot_data(result_05, "../../figures/chapter_NN/result_05.pdf")

set.seed(100)
result_06 <- lapply(1:40, function(x) learnNN_02(DT.model, 20, 'adam'))
saveRDS(result_06, file = "result_06.rds")
analyse_06 <- analyse_nn(result_06)
#plot_data(result_06, "../../figures/chapter_NN/result_06.pdf")

set.seed(100)
result_07 <- lapply(1:40, function(x) learnNN_02(DT.model, 37, 'adam'))
saveRDS(result_07, file = "result_07.rds")
analyse_07 <- analyse_nn(result_07)
#plot_data(result_07, "../../figures/chapter_NN/result_07.pdf")

set.seed(100)
result_08 <- lapply(1:40, function(x) learnNN_02(DT.model, 50, 'adam'))
saveRDS(result_08, file = "result_08.rds")
analyse_08 <- analyse_nn(result_08)
#plot_data(result_08, "../../figures/chapter_NN/result_08.pdf")


###################################################
############### Three hidden layers ############### 
###################################################

set.seed(100)
result_01_sgd <- lapply(1:40, function(x) learnNN_01(DT.model, 8, 'sgd'))
saveRDS(result_01_sgd, file = "result_01_sgd.rds")
#plot_data(result_01, "result_01.pdf")

set.seed(100)
result_02_sgd <- lapply(1:40, function(x) learnNN_01(DT.model, 20, 'sgd'))
saveRDS(result_02_sgd, file = "result_02_sgd.rds")
#plot_data(result_02, "result_02.pdf")

set.seed(100)
result_03_sgd <- lapply(1:40, function(x) learnNN_01(DT.model, 37, 'sgd'))
saveRDS(result_03_sgd, file = "result_03_sgd.rds")
#plot_data(result_03, "result_03.pdf")

set.seed(100)
result_04_sgd <- lapply(1:40, function(x) learnNN_01(DT.model, 50, 'sgd'))
saveRDS(result_04_sgd, file = "result_04_sgd.rds")
#plot_data(result_04, "result_04.pdf")

#################################################
############### Two hidden layers ###############
#################################################

set.seed(100)
result_05_sgd <- lapply(1:40, function(x) learnNN_02(DT.model, 8, 'sgd'))
saveRDS(result_05_sgd, file = "result_05_sgd.rds")
#plot_data(result_05, "result_05.pdf")

set.seed(100)
result_06_sgd <- lapply(1:40, function(x) learnNN_02(DT.model, 20, 'sgd'))
saveRDS(result_06_sgd, file = "result_06_sgd.rds")
#plot_data(result_06, "result_06.pdf")

set.seed(100)
result_07_sgd <- lapply(1:40, function(x) learnNN_02(DT.model, 37, 'sgd'))
saveRDS(result_07_sgd, file = "result_07_sgd.rds")
#plot_data(result_07, "result_07.pdf")

set.seed(100)
result_08_sgd <- lapply(1:40, function(x) learnNN_02(DT.model, 50, 'sgd'))
saveRDS(result_08_sgd, file = "result_08_sgd.rds")
#plot_data(result_08, "result_08.pdf")
