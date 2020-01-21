
############## Analyse the list of results from the neural network  ##############

analyse_nn <- function(results) {
  
  # bind all results from the runs into one long table together
  res_long <- rbindlist(results)
  # calculate the absoulte difference between predicted and true for all entries
  res_long[, abs_diff := abs(predicted - true)]
  # calculate the mean per year for the predicted values and bind together with true values
  res_agg <- res_long[, .(predicted = mean(predicted)), by = .(Year)] %>% cbind(results[[1]][, .(true)])
  # calculate relative and absolute differences per year
  res_agg[, abs_diff := abs(predicted - true) ]
  res_agg[, rel_diff := abs_diff / true ]
  
  Y25_mean_smoothed <- mean(res_agg[1:25, rel_diff], na.rm = TRUE)
  Y60_mean_smoothed <- mean(res_agg[1:61, rel_diff], na.rm = TRUE)
  
  Y25_mean <- res_long[Year %in% 1:25, .(sum(abs_diff, na.rm = TRUE) / sum(true, na.rm = TRUE))]
  Y60_mean <- res_long[, .(sum(abs_diff, na.rm = TRUE) / sum(true, na.rm = TRUE))]
  
  return(list(overview = res_agg,
              Y25_mean_smoothed = Y25_mean_smoothed,
              Y25_mean_smoothed = Y25_mean_smoothed,
              Y25_mean = Y25_mean,
              Y60_mean = Y60_mean))

}
