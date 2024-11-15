library(xtable)


# simulations parameters to change
source("estimators_mis.R")
NAME_OF_SIM <- "setting_1_mis"
SLIDER <- NULL
f <- f_optimal
# specify which estimators to use in simulation
# ordering of variables associated functions should match
var_order <- list("ipw_b", "ipw_nb", "gc_b", "gc_nb")
fun_order <- list(
  calculate_ipw_b, calculate_ipw_nb,
  calculate_gcomp_b, calculate_gcomp_nb
)

# simulation parameter that should not be changed
n_iter <- 2000 # number of simulation interrations, should be 2000
n_individuals <- c(2000, 500, 200) # number of individulas 2000,500,200
SEED <- 24091996
CORES_NOT_USED <- 10 # 10 for cluster, 1 for local machine



# help functions
est_clinical_util <- function(df, f, estimator) {
  return(mean(df$y) - estimator(df, f))
}

is_in_ci <- function(df, f, estimator, true_value, n_sample = 500) {
  n_row <- nrow(df)
  values <- rep(0, n_sample)
  for (i in 1:n_sample) {
    values[i] <- est_clinical_util(
      df = df[sample(n_row, n_row, replace = TRUE), ], # nolint: line_length_linter.
      f = f,
      estimator = estimator
    )
  }
  quantiles <- quantile(values, probs = c(0.025, 0.975))
  q_1 <- unname(quantiles[1])
  q_2 <- unname(quantiles[2])

  return(true_value > q_1 & true_value < q_2)
}

one_sim_step <- function(coverage = TRUE) {
  res <- generate_data_nb(n, slider = SLIDER) # nolint
  df <- res[[1]]
  k <- length(var_order)

  # return vector is structured in the following way, k number of estimators
  # k x est_clinical_utility, k x est_clinical_utility, index, mean, counterfactual, true_mean, true_counterfactual, true_clinical_utility # nolint
  return_vec <- rep(0, 2 * k + 6)

  return_vec[2 * k + 1] <- i
  return_vec[2 * k + 2] <- res[[2]]
  return_vec[2 * k + 3] <- res[[3]]
  return_vec[2 * k + 4] <- std_of_care
  return_vec[2 * k + 5] <- outcome_c
  return_vec[2 * k + 6] <- clinical_util

  for (j in seq_along(var_order)) {
    return_vec[j] <- est_clinical_util(
      df = df,
      f = f,
      estimator = fun_order[[j]]
    ) # nolint
  }
  if (coverage) {
    for (j in seq_along(var_order)) {
      return_vec[j + length(var_order)] <- is_in_ci(
        df = df,
        f = f,
        estimator = fun_order[[j]],
        true_value = clinical_util
      ) # nolint: line_length_linter.
    }
  }
  return(return_vec)
}

# -------run simulation---------

set.seed(SEED)

# calculate true value
big_df <- generate_data_nb(n = 10000000, slider = SLIDER)
std_of_care <- big_df[[2]]
outcome_c <- big_df[[3]]

clinical_util <- std_of_care - outcome_c

if (!file.exists("sim_data")) {
  dir.create("sim_data")
}

result_matrix <- as.matrix(var_order)

# _____parallel computing_____

library(foreach)
library(doParallel)

for (n in n_individuals) {
  cores <- detectCores()
  cl <- makeCluster(cores[1] - CORES_NOT_USED)
  registerDoParallel(cl)



  rm(final_matrix)
  final_matrix <- foreach(
    i = 1:n_iter,
    .combine = cbind,
    .packages = c("nnet")
  ) %dopar% {
    temp_matrix <- one_sim_step()

    temp_matrix
  }

  # stop cluster
  stopCluster(cl)

  # write simulated data into file for laster analysis
  df <- data.frame(final_matrix)
  rownames(df) <- c(
    paste0(var_order, "_estimand"),
    paste0(var_order, "_is_ci_correct"),
    "index",
    "mean",
    "counterfactual",
    "true_mean",
    "true_counterfactual",
    "true_clinical_util"
  )

  saveRDS(df, file = paste0("sim_data/", NAME_OF_SIM, "_", n, ".rds"))

  # calculate statistics of interest

  cat("Date:", Sys.time(), "\n")
  cat("\n Simulation with n = ", n, "\n") # nolint
  temp_result <- numeric(length(var_order))
  # Bias
  for (j in seq_along(var_order)) {
    temp_result[j] <- round(mean(final_matrix[j, ] - clinical_util) * 100, 3)
    cat("Bias * 100 of", var_order[[j]], ":", temp_result[j], "\n") # nolint
  }
  result_matrix <- cbind(result_matrix, temp_result)

  # Emperical Standard Deviation
  for (j in seq_along(var_order)) {
    temp_result[j] <- round(sd(final_matrix[j, ]) * 10, 3)
    cat("ESD * 10 of", var_order[[j]], ":", temp_result[j], "\n") # nolint
  }
  result_matrix <- cbind(result_matrix, temp_result)


  # coverage
  for (j in seq_along(var_order)) {
    temp_result[j] <- round(mean(final_matrix[length(var_order) + j, ]), 3)
    cat("Coverage of", var_order[[j]], ":", temp_result[j], "\n") # nolint
  }
  result_matrix <- cbind(result_matrix, temp_result)
}

cat("\n same information in partial latex syntax: \n")
print(xtable(result_matrix, type = "latex"))
