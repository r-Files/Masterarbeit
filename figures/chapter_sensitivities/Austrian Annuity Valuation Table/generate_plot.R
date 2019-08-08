############## load/install all required packages ##############
suppressWarnings(library("pacman"))
pacman::p_load(data.table, dplyr, rstudioapi, ggplot2, scales)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read the data starting from line with "Alter"
annTable <- fread("QxAV05RE_Unisex.txt",
                  skip = "Alter")

annTable %>% setnames(old = c("Alter", "q(2001)Einzel(u)"),
                      new = c("age", "qx"))

p2 <- ggplot(annTable, aes(x = age, y = qx)) + geom_line() +
  scale_x_continuous(limits = c(0, 120)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title   = "Austrian Annuity Valuation Table AVÖ 2005R Unisex",
       caption = "Source: Actuarial Association of Austria",
       y       = expression("log "*"q"[x]*"(2001)"))
  
# save the plot as pdf
ggsave("avoe_2005R_unisex.pdf")
