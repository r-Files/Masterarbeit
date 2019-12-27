library(ggplot2)

fun_sigmoid <- function(x) 1 / (1 + exp(-x))
fun_tanh <- function(x) (exp(x) - exp(-x)) / (exp(x) + exp(-x))
fun_relu <- function(x) pmax(0,x)
fun_softplus <- function(x) log(1 + exp(x))


base <- ggplot(data.frame(x = c(-3, 3)), aes(x))

base + stat_function(fun = fun_sigmoid, aes(colour = "Sigmoid")) +
  stat_function(fun = fun_tanh, aes(colour = "tanh")) +
  stat_function(fun = fun_softmax, aes(colour = "Softplus")) +
  stat_function(fun = fun_relu, aes(colour = "ReLu")) +
  scale_colour_manual("Activation functions", values = c("red", "blue", "green", "orange"))
