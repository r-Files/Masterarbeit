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


DT.input <- prepare_input()
DT.ouput <- prepare_output(DT.input$DGO_SPCODE_I)
DT.model <- split_data(DT.input, DT.ouput, 0.80)

############### Three hidden layers ############### 

result_01 <- learnNN_01(DT.model)
ggplot(data = melt(result_01, id = "Year", value.name = "reserve") ,
       aes(x=Year, y = reserve, colour = variable)) +
  geom_line()

result_02 <- learnNN_02(DT.model)
ggplot(data = melt(result_02, id = "Year", value.name = "reserve") ,
       aes(x=Year, y = reserve, colour = variable)) +
  geom_line()

result_03 <- learnNN_03(DT.model)
ggplot(data = melt(result_02, id = "Year", value.name = "reserve") ,
       aes(x=Year, y = reserve, colour = variable)) +
  geom_line()

result_04 <- learnNN_04(DT.model)
ggplot(data = melt(result_04, id = "Year", value.name = "reserve") ,
       aes(x=Year, y = reserve, colour = variable)) +
  geom_line()

############### Two hidden layers ###############

result_05 <- learnNN_05(DT.model)
ggplot(data = melt(result_05, id = "Year", value.name = "reserve") ,
       aes(x=Year, y = reserve, colour = variable)) +
  geom_line()

result_06 <- learnNN_06(DT.model)
ggplot(data = melt(result_06, id = "Year", value.name = "reserve") ,
       aes(x=Year, y = reserve, colour = variable)) +
  geom_line()

result_07 <- learnNN_07(DT.model)
ggplot(data = melt(result_07, id = "Year", value.name = "reserve") ,
       aes(x=Year, y = reserve, colour = variable)) +
  geom_line()

result_08 <- learnNN_08(DT.model)
ggplot(data = melt(result_08, id = "Year", value.name = "reserve") ,
       aes(x=Year, y = reserve, colour = variable)) +
  geom_line()


test_data_long[, Year := as.factor(Year)]



