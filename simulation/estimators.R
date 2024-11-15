library(nnet)


f_optimal <- function(df) {
  t <- numeric(nrow(df))

  t[df$z_1 < 0.3] <- 1
  t[df$z_1 > 0.7] <- 3
  t[(df$z_1 >= 0.3) & (df$z_1 <= 0.7)] <- 2

  return(t)
}

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

generate_data_nb <- function(n = 1000, print_summary = FALSE, slider = NULL) {
  z_1 <- sample(seq(from = 0, to = 1, by = 0.1), n, replace = TRUE)
  z_2 <- rbinom(n, size = 1, prob = 0.7)

  p_1 <- numeric(n)
  p_2 <- numeric(n)
  p_3 <- numeric(n)
  t <- numeric(n)

  # standard of care distibution via
  if (is.null(slider)) {
    p_1 <- plogis(-z_1)
    p_2 <- plogis(-2 + z_1 + z_2)
    p_3 <- 1 - p_1 - p_2

    t <- mapply(
      function(x, y, z) sample(c(1, 2, 3), size = 1, replace = TRUE, prob = c(x, y, z)), # nolint: line_length_linter.
      p_1, p_2, p_3
    )
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

    t <- mapply(
      function(x, y, z) sample(c(1, 2, 3), size = 1, replace = TRUE, prob = c(x, y, z)), # nolint: line_length_linter.
      p_1, p_2, p_3
    )
  }
  prob_default <- p_default(t, z_1, z_2)
  y <- rbinom(n, prob = prob_default, size = 1)

  std_care_avg <- mean(prob_default)

  t_c <- f(data.frame(z_1))
  prob_default_c <- p_default(t_c, z_1, z_2)
  counterfactual_avg <- mean(prob_default_c)

  if (print_summary) {
    cat("average outcome:", std_care_avg, "\n")
    cat("counterfactual outcome:", counterfactual_avg, "\n")
  }
  sim_data <- data.frame(z_1, z_2, t, y, t_c) # nolint: line_length_linter.
  return(list(sim_data, std_care_avg, counterfactual_avg))
}


calculate_ipw_b <- function(df, f) {
  w <- f(df) == df$t
  df$z_1 <- as.factor(df$z_1)

  model <-
    glm(w ~ z_1 * z_2, family = binomial(link = "logit"), data = df) # nolint
  weights <-
    predict(model, df[w, c("z_1", "z_2")], type = "response")
  expected_y <- sum(df$y[w] / weights) / nrow(df)
  return(expected_y)
}

calculate_ipw_nb <- function(df, f) {
  w <- f(df) == df$t
  df$z_1 <- as.factor(df$z_1)


  model <- multinom(t ~ z_1 * z_2, data = df, trace = FALSE) # nolint: line_length_linter.
  all_weights <-
    predict(model, newdata = df[w, c("z_1", "z_2")], type = "probs") # nolint


  # select the weights the treatment would choose
  weights <- numeric(nrow(all_weights))
  for (i in seq_len(nrow(all_weights))) {
    weights[i] <- all_weights[i, df$t[w][i]]
  }

  expected_y <- sum(df$y[w] / weights) / nrow(df)
  expected_y
  return(expected_y)
}


calculate_gcomp_b <- function(df, f) {
  w <- f(df) == df$t
  df$z_1 <- as.factor(df$z_1)
  df$proposed_t <- w

  model <-
    glm(df$y ~ z_1 * z_2 * proposed_t, # nolint: line_length_linter.
      family = binomial(link = "logit"),
      data = df
    )

  df$proposed_t <- TRUE

  weights <-
    predict(model, df[c("z_1", "z_2", "proposed_t")], type = "response") # nolint: line_length_linter.
  expected_y <- sum(weights) / nrow(df)
  return(expected_y)
}


calculate_gcomp_nb <- function(df, f) {
  df$t1 <- as.integer(df$t == 1)
  df$t2 <- as.integer(df$t == 2)
  df$t3 <- as.integer(df$t == 3)

  model <-
    glm(df$y ~ 0 + t1 + t1:(z_1 + z_2) + t2 + t2:(z_2) + t3 + t3:(z_1 + z_2), # nolint: line_length_linter.
      data = df
    )

  df$t1 <- as.integer(f(df) == 1)
  df$t2 <- as.integer(f(df) == 2)
  df$t3 <- as.integer(f(df) == 3)

  weights <-
    predict(model, newdata = df, type = "response") # nolint: line_length_linter.
  expected_y <- sum(weights) / nrow(df)
  return(expected_y)
}

calculate_dr <- function(df, f) {
  # implcit assumtion: g_comp_nb estiamtes with non categorical z_1 but ipw_b does # nolint

  # estimating outcomes via g_comp_nb
  df$t1 <- as.integer(df$t == 1)
  df$t2 <- as.integer(df$t == 2)
  df$t3 <- as.integer(df$t == 3)


  model <-
    glm(df$y ~ 0 + t1 + t1:(z_1 + z_2) + t2 + t2:(z_2) + t3 + t3:(z_1 + z_2), # nolint: line_length_linter.
      data = df
    )

  df$t1 <- as.integer(f(df) == 1)
  df$t2 <- as.integer(f(df) == 2)
  df$t3 <- as.integer(f(df) == 3)

  y_hat <-
    predict(model, df[c("z_1", "z_2", "t1", "t2", "t3")], type = "response") # nolint: line_length_linter.

  # estemating propensity scores via ipw_b
  w <- f(df) == df$t
  df$z_1 <- as.factor(df$z_1)

  model <-
    glm(w ~ z_1 * z_2, family = binomial(link = "logit"), data = df) # nolint
  propensity_score <-
    predict(model, df[w, c("z_1", "z_2")], type = "response")

  # calculation of double robust estimator as defined in paper
  expected_y <- sum(y_hat) / nrow(df) + sum((1 / propensity_score) * (df$y[w] - y_hat[w])) / nrow(df) # nolint: line_length_linter.

  return(expected_y)
}
