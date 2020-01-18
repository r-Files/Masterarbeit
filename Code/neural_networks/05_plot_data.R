
############## Plot the data from the neural network  ##############

plot_data <- function(input, filename) {
  
  # melt the input data from wide to long format
  data_long <- melt(input, id = "Year", value.name = "reserve")
  
  # generate a suitable plot
  ggplot(data = data_long ,
         aes(x=Year, y = reserve, colour = variable)) +
    geom_line(size = 0.75) +
    theme(legend.title    = element_blank(),
          legend.position = "bottom")
  
  # save the plot as .pdf
  ggsave(filename, 
         plot = last_plot(), # or give ggplot object name as in myPlot,
         width = 11.69, height = 8.27, 
         units = "in", # other options c("in", "cm", "mm"), 
         dpi = 300)

}
