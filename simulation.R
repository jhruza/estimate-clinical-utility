library(ggplot2)
source("estimator_2d.R")


#' override decision function, for documentation see estimator_2d.R
decision_function <- function(df) {
  t <- df$z_2 + df$z_1 > 1
  return(t)
}

# simulation paramters
n_iter <- 1000 #number of simulation interrations
n <- 1000 #number of individulas per simulation

# empty lists to be filled with results during summulation
res_truth <- list()
res_naive_ipw <- list()
res_ipw <- list()
res_ipw_exact <- list()
res_g_c <- list()
res_naive_g_c <- list()
res_dr <- list()
res_dr_exact <- list()


for (i in 1:n_iter) {

  res <- generate_data(n = n, f = decision_function)
  df <- res[[1]]
  res_truth <- c(res_truth, res[[3]])
  res_naive_ipw <- c(res_naive_ipw, calculate_naive_ipw(df, decision_function))
  res_ipw <- c(res_ipw, calculate_ipw(df, f = decision_function, estimated_weights = TRUE)) # nolint: line_length_linter.
  res_ipw_exact <- c(res_ipw_exact, calculate_ipw(df, f = decision_function, estimated_weights = FALSE)) # nolint: line_length_linter.
  res_g_c <- c(res_g_c, calculate_g_comp(df, decision_function, estimated_weights = FALSE)) # nolint: line_length_linter.
  res_naive_g_c <- c(res_naive_g_c, calculate_naive_g_comp(df, decision_function)) # nolint: line_length_linter.
  res_dr <- c(res_dr, calculate_dr(df, decision_function)) # nolint: line_length_linter.
  res_dr_exact <- c(res_dr_exact, calculate_dr(df, decision_function, estimated_weights = FALSE)) # nolint: line_length_linter.

  #print progress
  if (i %% (n_iter / 100) == 0) {
    cat("\r", i / n_iter * 100, "%")
    flush.console()
  }
}
cat("\r")
flush.console()

# plot
# Combine data into a data frame for a plot
data <- data.frame(
  value = c(unlist(res_truth), unlist(res_naive_ipw), unlist(res_ipw), unlist(res_ipw_exact), unlist(res_g_c), unlist(res_naive_g_c), unlist(res_dr), unlist(res_dr_exact)), # nolint: line_length_linter.
  group = rep(c("truth",
                "naive IPW",
                "IPW",
                "IPW exact",
                "g computation",
                "naive g copmutation",
                "double robust",
                "double robust exact"),
              each = n_iter)
)

mean_values <- aggregate(value ~ group, data, mean)

# Create a ggplot with shared axes and mean lines
ggplot(data, aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.7, bins = 100) +
  geom_vline(data = mean_values, aes(xintercept = value, color = group),
             linetype = "dashed", size = 1) +
  facet_grid(rows = vars(group), scales = "free") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1))

mean_values


ggplot() +
  geom_line(aes(x =  1:n_iter,
                y = cumsum(unlist(res_ipw)) / 1:n_iter),
            color = "red") +
  geom_line(aes(x = 1:n_iter,
                y = cumsum(unlist(res_ipw_exact)) / 1:n_iter),
            color = "blue") +
  geom_line(aes(x = 1:n_iter,
                y = cumsum(unlist(res_truth)) / 1:n_iter),
            color = "black") +
  geom_line(aes(x = 1:n_iter,
                y = cumsum(unlist(res_naive_ipw)) / 1:n_iter),
            color = "#009600") +
  geom_line(aes(x = 1:n_iter,
                y = cumsum(unlist(res_dr_exact)) / 1:n_iter),
            color = "#e108ab") +
  geom_line(aes(x = 1:n_iter,
                y = cumsum(unlist(res_g_c)) / 1:n_iter),
            color = "#08e1da") +
  geom_line(aes(x = 1:n_iter,
                y = cumsum(unlist(res_dr)) / 1:n_iter),
            color = "#e16208") +
  ylab("Values") + xlab("date")

# single run
df <- generate_data(n = 10000, f = decision_function)[[1]] # nolint: line_length_linter.
mean(df$y_c)
mean(df$y)
mean(df$y[decision_function(df) == df$treat])

plot(z_1 ~ z_2,
     xlab = "z2",
     ylab = "z1",
     main = "weights of used datapoints",
     col = ifelse(df$y_c, "red", "green"),
     pch = ifelse(df$treat_c, 3, 1),
     data = df)

calculate_naive_ipw(df, decision_function, plot = TRUE)
calculate_naive_ipw(df, decision_function, plot = TRUE, estimated_weights = FALSE) # nolint: line_length_linter.
calculate_ipw(df, decision_function, plot =  TRUE, estimated_weights = TRUE)
calculate_ipw(df, decision_function, plot =  TRUE, estimated_weights = FALSE)
calculate_g_comp(df, decision_function, estimated_weights = TRUE)
calculate_g_comp(df, decision_function, estimated_weights = FALSE)
calculate_naive_g_comp(df, decision_function)
calculate_dr(df, decision_function, estimated_weights = FALSE)
