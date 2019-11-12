suppressWarnings(library("pacman"))
pacman::p_load(data.table,
               ggplot2,
               dplyr,
               nnls,
               matrixStats,
               microbenchmark,
               Rcpp,
               pracma, # for Householder
               Rfast,
               matlib) #bit64

Hauptverzeichnis <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/")
setwd (Hauptverzeichnis)

#
# Test accuracy of QR-decomposition
#

method_one <- function(A, b) {
  solve((t(A) %*% A)) %*% t(A) %*% b
}

method_two <- function(A, b) {
  solve(t(A) %*% A, t(A) %*% b)
}

method_three <- function(A, b) {
  solve(crossprod(A), crossprod(A,b))
}

method_four <- function(A, b) {
  deco <- qr(A, LAPACK = FALSE)
  solve(qr.R(deco)) %*% t(qr.Q(deco)) %*% b
}

method_five <- function(A, b) {
  deco <- qr(A, LAPACK = FALSE)
  solve(qr.R(deco), t(qr.Q(deco)) %*% b)
}

method_six <- function(A, b) {
  deco <- qr(A, LAPACK = FALSE)
  backsolve(qr.R(deco), t(qr.Q(deco)) %*% b)
}

method_seven <- function(A, b) {
  qr.solve(qr(A, LAPACK = FALSE), b)
}






# make the results reproducable
set.seed(100)

m <- 5000
# generate Ax = b
# start with A being mxm and solve least square problem for
# m x 50, m x 100, m x 150, m x 200, m x 250, ...
A <- matrix(runif(m * m, 0, 100000), nrow = m)
x <- runif(m, 0, 1000)

#cols <- lapply(c(10, 50, 100, 150, 200, 300, 400, 500, 750, 1000), seq_len)
number_cols <- c(10, 50, 100, 150, 200, 300, 400, 500, 750, 1000, 1500, 2000, 2500, 3000)
cols <- lapply(number_cols, seq_len)

difference <- lapply(cols, function(c) {
  A_P <- A[, c]
  x_P <- x[c]
  b_P <- A_P %*% x_P
  
  meth_one <- method_one(A_P, b_P)
  meth_two <- method_two(A_P, b_P)
  meth_thr <- method_three(A_P, b_P)
  meth_fou <- method_four(A_P, b_P)
  meth_fiv <- method_five(A_P, b_P)
  meth_six <- method_six(A_P, b_P)
  meth_sev <- method_seven(A_P, b_P)
  
  errors <- list(
    "method_one"   = abs(x_P - meth_one) %>% mean(),
    "method_two"   = abs(x_P - meth_two) %>% mean(),
    "method_three" = abs(x_P - meth_thr) %>% mean(),
    "method_four"  = abs(x_P - meth_fou) %>% mean(),
    "method_five"  = abs(x_P - meth_fiv) %>% mean(),
    "method_six"   = abs(x_P - meth_six) %>% mean(),
    "method_seven" = abs(x_P - meth_sev) %>% mean()
  )
  
  residuals <- list(
    "method_one"   = abs(A_P %*% meth_one - b_P) %>% sum(),
    "method_two"   = abs(A_P %*% meth_two - b_P) %>% sum(),
    "method_three" = abs(A_P %*% meth_thr - b_P) %>% sum(),
    "method_four"  = abs(A_P %*% meth_fou - b_P) %>% sum(),
    "method_five"  = abs(A_P %*% meth_fiv - b_P) %>% sum(),
    "method_six"   = abs(A_P %*% meth_six - b_P) %>% sum(),
    "method_seven" = abs(A_P %*% meth_sev - b_P) %>% sum()
  )
  
  list("err" = errors,
       "res" = residuals)
})

names(difference) <- number_cols



micro <- lapply(cols, function(c) {
  A_P <- A[, c]
  x_P <- x[c]
  b_P <- A_P %*% x_P
  
  timing <- microbenchmark(
    times = 10,
    method_one(A_P, b_P) ,
    method_two(A_P, b_P),
    method_three(A_P, b_P),
    method_four(A_P, b_P),
    method_five(A_P, b_P),
    method_six(A_P, b_P),
    method_seven(A_P, b_P)
  )
  
  tmp <- data.table(time = timing$time,
                    form = timing$expr)
  tmp[, .(avg = mean(time) / 1e9,
          med = median(time) / 1e9),
      by = .(form)][order(form)]
})

names(micro) <- number_cols

micro <- rbindlist(micro, idcol = "col_cnt")
micro$col_cnt <- as.numeric(micro$col_cnt)

levels(micro$form)[levels(micro$form) == "method_one(A_P, b_P)"] <- "method one"
levels(micro$form)[levels(micro$form) == "method_two(A_P, b_P)"] <- "method two"
levels(micro$form)[levels(micro$form) == "method_three(A_P, b_P)"] <- "method three"
levels(micro$form)[levels(micro$form) == "method_four(A_P, b_P)"] <- "method four"
levels(micro$form)[levels(micro$form) == "method_five(A_P, b_P)"] <- "method five"
levels(micro$form)[levels(micro$form) == "method_six(A_P, b_P)"] <- "method six"
levels(micro$form)[levels(micro$form) == "method_seven(A_P, b_P)"] <- "method seven"

# define a column to distinguish between methods which use QR and not
micro[, QR := "QR used"]
# methods one to three don't use QR decomposition
micro[form %in% c("method one", "method two", "method three"), QR := "QR not used"]


p_runtime <-
  ggplot(micro,
         aes(
           x = col_cnt,
           y = avg,
           group = form,
           shape = form,
           color = form
         )) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values = 1:7) +
  facet_grid(. ~ QR) +
  labs(x = "Number of columns",
       y = "Average computation time in seconds",
       title = "Computational time least squares problem") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
p_runtime


saveRDS(micro, "micro.rds")
saveRDS(difference, "difference.rds")


micro <- readRDS("micro.rds")
difference <- readRDS("difference.rds")

###############
## errors
###############

# extract the errors from the list and melt it into a long format
errors_for_plotting <- sapply(difference, function(x) {
  x$err %>% unlist()
}) %>% as.data.table(keep.rownames = "method") %>%
  melt(id.vars = "method",
       variable.name = "col_cnt",
       value.name = "error")

# make sure that the levels are ordered in a way we want it to
errors_for_plotting[, method := factor(
  method,
  levels = c(
    "method_one",
    "method_two",
    "method_three",
    "method_four",
    "method_five",
    "method_six",
    "method_seven"
  )
)]

# plot the errors for each number of columns
p_error <-
  ggplot(errors_for_plotting[!(col_cnt %like% "^(10|50)"), ], # kick out what's like 10.. or 50..
         aes(x = method,y = error)) +
  geom_point() +
  facet_wrap(~ col_cnt, scales = "free_y") +
  labs(title = "Average errors per column count and method") +
  ylab(expression(paste("Average absolute error between ", s[true], " and ", s[LS]))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())
p_error


###############
## residuals
###############

residuals_for_plotting <- sapply(difference, function(x) {
  x$res %>% unlist()
}) %>% as.data.table(keep.rownames = "method") %>%
  melt(id.vars = "method",
       variable.name = "col_cnt",
       value.name = "residuals")

# make sure that the levels are ordered in a way we want it to
residuals_for_plotting[, method := factor(
  method,
  levels = c(
    "method_one",
    "method_two",
    "method_three",
    "method_four",
    "method_five",
    "method_six",
    "method_seven"
  )
)]

p_resid <-
  ggplot(residuals_for_plotting[!(col_cnt %like% "^(10|50)"), ], # kick out what's like 10.. or 50..
         aes(x = method,y = residuals)) +
  geom_point() +
  facet_wrap(~ col_cnt, scales = "free_y") +
  labs(title = "Sum of residuals per column count and method") +
  ylab(expression(paste("Sum of absolute residuals: ", A * s[LS] - b))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())
p_resid
