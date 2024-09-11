#' calculate_ipw_b
#'
#' Calculates the average outcome had the given dataset been treated with the given treatment regime.
#'
#' @param df The dataset of observed data, must contain columns treatment: "t", outcome: "y" and covariates. (data.frame)
#' @param f The treatment regime. (function)
#' @param covariates Optional. The covariates used for propensity score. Must be a subset of covariates in df. (list)
#' @param propensity_model_equation Optional. The RHS equation for the propensity model, of the form "covariate_1 + covariate_2" (string)
#'
#' @return The calculated inverse probability weights.
#'
#' @export
calculate_ipw_b <- function(df, f, covariates = NULL, propensity_model_equation = NULL) {
  stopifnot(is.data.frame(df), "t" %in% names(df), "y" %in% names(df))

  # define covariates used for the propensity score and generate formula for propensity score model
  if (is.null(propensity_model_equation)) {
    covariates <- setdiff(names(df), c("t", "y"))
    formula_string <- paste("w ~", paste(covariates, collapse = " * "))
  } else {
    # check if propensity_model_equation is a string
    stopifnot(class(propensity_model_equation) == "string")
    # check that covariates are in df
    stopifnot(all(covariates %in% names(df)))

    formula_string <- paste("w ~", propensity_model_equation)
  }

  formula <- as.formula(formula_string)
  w <- f(df) == df$t

  model <-
    glm(formula, family = binomial(link = "logit"), data = df)
  weights <-
    predict(model, df[w, covariates], type = "response")
  expected_y <- sum(df$y[w] / weights) / nrow(df)

  return(list(outcome = expected_y, prop_model = model))
}
#' Calculate Clinical Utility
#'
#' This function calculates the clinical utility of two treatment assignment functions `f` and `g`
#' using inverse probability weighting (IPW). Optionally, it can also estimate the variance of the
#' clinical utility difference.
#'
#' @param df A data frame containing the dataset. It must include columns "t" (treatment) and "y" (outcome).
#' @param f A treatement regime that assigns treatment based on the covariates in `df`.
#' @param g Another treatment regime that assigns treatment based on the covariates in `df`.
#' @param est_variance A logical value indicating whether to estimate the variance of the clinical utility
#' difference. Default is TRUE.
#'
#' @return A list containing:
#' \item{clinical_util}{The difference in clinical utility between the two treatment assignment functions `f` and `g`.}
#' \item{variance}{(Optional) The estimated variance of the clinical utility difference, if `est_variance` is TRUE.}
#'
#' @export
clinical_utility <- function(df, f, g, est_variance = TRUE) {
  n <- nrow(df)
  ipw_f <- calculate_ipw_b(df, f)
  ipw_g <- calculate_ipw_b(df, g)
  clinical_util <- ipw_f$outcome - ipw_g$outcome
  if (est_variance) {
    covariates <- df[, setdiff(names(df), c("t", "y"))]

    # indicators if assigned treatments follows observed treatment
    w_f <- f(df) == df$t
    w_g <- g(df) == df$t

    # propensity scores for f
    ps_f <- predict(ipw_f$prop_model, covariates, type = "response")
    # propensity scores for g
    ps_g <- predict(ipw_g$prop_model, covariates, type = "response")

    # assume model matrix is the same for ipw_f and ipw_g
    model_mat <- model.matrix(ipw_f$prop_model)
    n_features <- ncol(model_mat)
    A <- colSums(model_mat * w_f * df$y * (1 - (1 / ps_f))) / n
    B <- colSums(model_mat * w_g * df$y * (1 - (1 / ps_g))) / n

    partial_M2 <- t(model_mat) %*% (model_mat * (-ps_f + ps_f^2))
    partial_M3 <- t(model_mat) %*% (model_mat * (-ps_g + ps_g^2))

    partial_M2 <- partial_M2 / n
    partial_M3 <- partial_M3 / n

    partial_M2_inv <- solve(partial_M2)
    partial_M3_inv <- solve(partial_M3)

    # create bread matrix for vegetarian sandwich estimator
    bread_matrix <- rbind(
      cbind(-1, A %*% partial_M2_inv, B %*% partial_M3_inv),
      cbind(0, partial_M2_inv, matrix(0, nrow = n_features, ncol = n_features)),
      cbind(0, matrix(0, nrow = n_features, ncol = n_features), partial_M3_inv)
    )

    # create cheese matrix for vegetarian sandwich estimator
    cheese_matrix <- matrix(0, nrow = 2 * n_features + 1, ncol = 2 * n_features + 1)
    for (i in 1:n) {
      c_i <- as.numeric(model_mat[i, ])
      M1 <- as.numeric(w_f[i] * df$y[i] / ps_f[i] - w_g[i] * df$y[i] / ps_g[i] - clinical_util)
      M2 <- c_i * (w_f[i] - ps_f[i])
      M3 <- c_i * (w_g[i] - ps_g[i])
      m_i <- c(M1, M2, M3)
      cheese_matrix <- cheese_matrix + outer(m_i, m_i, "*")
    }
    cheese_matrix <- cheese_matrix / n
    cov_matrix <- bread_matrix %*% cheese_matrix %*% t(bread_matrix)
    variance <- cov_matrix[1, 1] / n
  } else {
    variance <- NULL
  }
  return(list(outcome = clinical_util, variance = variance))
}

#' Remove factors from a data frame
#'
#' This function takes a data frame as input and removes all factor columns from it.
#' It converts the factor columns to numeric columns by converting the factor levels to their corresponding numeric values.
#'
#' @param df A data frame.
#'
#' @return A modified data frame with factor columns converted to numeric columns.
rm_factors <- function(df) {
  factor_cols <- sapply(df, is.factor)
  df[factor_cols] <- lapply(df[factor_cols], function(x) as.numeric(as.character(x)))
  return(df)
}

#' Optimal Treatment Assignment
#'
#' This function assigns an optimal treatment for the data generated by generate_data_nb().
#'
#' @param df A data frame that contains the column `z_1`.
#' @return A numeric vector representing the assigned treatments.
f_optimal <- function(df) {
  df <- rm_factors(df)
  t <- numeric(nrow(df))
  t[df$z_1 < 0.3] <- 1
  t[df$z_1 > 0.7] <- 3
  t[(df$z_1 >= 0.3) & (df$z_1 <= 0.7)] <- 2

  return(t)
}

#' Default Probability Calculation
#'
#' This function calculates the default probabilities based on treatment and covariates.
#'
#' @param t A numeric vector representing the assigned treatments.
#' @param z_1 A numeric vector representing the first covariate.
#' @param z_2 A numeric vector representing the second covariate.
#' @return A numeric vector representing the calculated probabilities.
p_default <- function(t, z_1, z_2) {
  condition_1 <- (t == 1)
  condition_2 <- (t == 2)
  condition_3 <- (t == 3)

  p <- numeric(length(t))

  p[condition_1] <- 0.5 + 0.5 * z_1[condition_1] - 0.5 * z_2[condition_1]
  p[condition_2] <- 0.65 - 0.5 * z_2[condition_2]
  p[condition_3] <- 1 - 0.5 * z_1[condition_3] - 0.5 * z_2[condition_3]

  return(p)
}
#' Generate Data with Non-Binary Covariates
#'
#' This function generates a dataset with non-binary covariates and assigns treatments based on the covariates.
#' For more details see Appendix
#'
#' @param n An integer representing the number of samples to generate. Default is 10000.
#' @param f A treatment regime for which the true counteractual average outcome is calculated. Default is NULL.
#' @param print_summary A logical indicating whether to print a summary of the generated data. Default is FALSE.
#' @param slider A variable to adjust the propability to get the optimal treatment according to the optimal
#' decision function. slider ->infinity: prop of being treatet opimal=1.
#' slider-> -infinity prop of being treatet opimal=0. Default is NULL.
#'
#'  @return A list containing:
#' \item{sim_data}{A data frame containing the simulated population}
#' \item{std_care_avg}{A numeric value representing the average standard of care.}
#' \item{counterfactual_avg}{A numeric value representing the average counterfactual outcome according to f.}
#' @examples
#' generate_data_nb(n = 1000, print_summary = TRUE)
generate_data_nb <- function(n = 10000, f = NULL, print_summary = FALSE, slider = NULL) {
  z_1 <- sample(seq(from = 0, to = 1, by = 0.1), n, replace = TRUE)
  z_2 <- rbinom(n, size = 1, prob = 0.7)

  # set up empty vectors for the probability of three treatments
  p_1 <- numeric(n)
  p_2 <- numeric(n)
  p_3 <- numeric(n)
  t <- numeric(n)

  # define standard of care distibution
  if (is.null(slider)) {
    p_1 <- plogis(-z_1)
    p_2 <- plogis(-2 + z_1 + z_2)
    p_3 <- 1 - p_1 - p_2
  } else {
    idx_1 <- f_optimal(data.frame(z_1)) == 1
    idx_2 <- f_optimal(data.frame(z_1)) == 2
    idx_3 <- f_optimal(data.frame(z_1)) == 3

    p_1[!idx_1] <- plogis(-z_1[!idx_1] - slider)
    p_2[!idx_2] <- plogis(-2 + z_1[!idx_2] + z_2[!idx_2] - slider)
    p_3[!idx_3] <- plogis(2 * z_1[!idx_3] - z_2[!idx_3] - slider)

    p_1[idx_1] <- plogis(slider)
    p_2[idx_2] <- plogis(slider)
    p_3[idx_3] <- plogis(slider)
  }
  # set treatment according to distribution
  t <- mapply(
    function(x, y, z) sample(c(1, 2, 3), size = 1, replace = TRUE, prob = c(x, y, z)),
    p_1, p_2, p_3
  )


  prob_default <- p_default(t, z_1, z_2)
  y <- rbinom(n, prob = prob_default, size = 1)

  std_care_avg <- mean(prob_default)

  if (!is.null(f)) {
    t_c <- f(data.frame(z_1))
    prob_default_c <- p_default(t_c, z_1, z_2)
    counterfactual_avg <- mean(prob_default_c)
  } else {
    t_c <- NA
    counterfactual_avg <- NA
  }

  if (print_summary) {
    cat("average outcome:", std_care_avg, "\n")
    cat("counterfactual outcome:", counterfactual_avg, "\n")
  }
  sim_data <- data.frame(z_1 = as.factor(z_1), z_2, t, y, t_c)
  return(list(sim_data = sim_data, std_care_avg = std_care_avg, counterfactual_avg = counterfactual_avg))
}



# example calculations:

f_1 <- function(df) rep(1, nrow(df))
f_2 <- function(df) rep(2, nrow(df))
f_3 <- function(df) rep(3, nrow(df))


data <- generate_data_nb(f = f_2, n = 500, slider = -1, print_summary = TRUE)
data$std_care_avg
data$counterfactual_avg
df <- data$sim_data[, c("y", "t", "z_1", "z_2")]
a <- calculate_ipw_b(df, f_2)
a$outcome

clinical_utility(df, f_2, f_optimal)
clinical_utility(df, f_optimal, f_2)
n_loop <- 500
n <- 500
results <- numeric(n_loop)
variances <- numeric(n_loop)
for (i in 1:n_loop) {
  df <- generate_data_nb(f = f_optimal, n = n, slider = -10)$sim_data[, c("y", "t", "z_1", "z_2")]
  cu <- clinical_utility(df, f_optimal, f_1, est_variance = TRUE)
  results[i] <- cu$outcome
  variances[i] <- cu$variance
}

hist(variances)
abline(v = var(unlist(results)))
var(results)
mean(variances)
mean(variances) / var(results)

density_results <- density(results)
plot(density_results, main = "Density Plot of Results", xlab = "Outcome", col = "blue")
lines(density_results, col = "blue")
x <- seq(-1, 1, length.out = 1000)
y <- dnorm(x, mean = mean(results), sd = sqrt(mean(variances)))
lines(x, y, col = "red")



# random generated data example:
df <- data.frame(t = sample(0:1, 100, replace = TRUE), y = rnorm(100), x1 = rnorm(100), x2 = rnorm(100))
f <- function(df) {
  ifelse(df$x1 > 0, 1, 0)
}
g <- function(df) {
  ifelse(df$x2 > 0, 1, 0)
}



n_loop <- 5000
n <- 1000
results <- numeric(n_loop)
variances <- numeric(n_loop)
for (i in 1:n_loop) {
  df <- data.frame(t = sample(0:1, n, replace = TRUE), y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n))
  cu <- clinical_utility(df, f, g, est_variance = TRUE)
  results[i] <- cu$outcome
  variances[i] <- cu$variance
}

hist(variances)
abline(v = var(unlist(results)))
mean(variances) / var(results)

density_results <- density(results)
plot(density_results, main = "Density Plot of Results", xlab = "Outcome", col = "blue")
lines(density_results, col = "blue")
x <- seq(-1, 1, length.out = 1000)
y <- dnorm(x, mean = mean(results), sd = sqrt(mean(variances)))
lines(x, y, col = "red")
