
############## Read and prepare the output-data ##############

prepare_output <- function(rowID) {
  
  # read the raw-data output data
  DT.output_raw <- fread("RUN_986/DMWV__~ACT_W_IF_DMWV__~Base.gro",
                     sep = ",",
                     header = TRUE,
                     skip = "DGO_SPCODE_I")
  
  # select the unique policy-ID and MATH_RES_IF (acturial reserve) as output 
  DT.output_sel <- DT.output_raw[, .(DGO_SPCODE_I, MATH_RES_IF)]
  # free space 
  rm(DT.output_raw)
  # select only those rows corresponding to the unique policy-ID's used for the input
  DT.output_sel <- DT.output_sel[DGO_SPCODE_I %in% rowID, ]
  # add increasing year-column per DGO_SPCODE 
  DT.output_sel[, year := 0:(.N - 1), by = .(DGO_SPCODE_I)]
  # Create Table where the columns are the development of the reserve
  # and the rows are the different contracts.
  DT.output <- dcast(DT.output_sel, DGO_SPCODE_I ~ year , value.var = "MATH_RES_IF", fill = 0)
  
  return(DT.output)
}