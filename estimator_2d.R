#' Decision function that describes the intervention for binary treatments 1d
#'
#' @param df dataframe with column z
#' @return Vector of treatments f(z)
#' @examples
#' decision_function_binary(df)
decision_function <- function(df) {
  t <- df$z_2 + df$z_1 < 1
  return(t)
}

#' Calcultes expit function defined by expit(x) = 1/(1+exp(-x))
#' @param x Input value.
expit <- function(x) {
  return(1 / (1 + exp(-x)))
}

#' Calcultes expit function defined by expit(x) = 1/(1+exp(-x))
#' @param x Input value.
lossfunction_treat <- function(z_1, z_2) {
  #p1 = (z_1, z_2)
  #p2 = (z_2, z_1)
  # vector from p1 to p2
  vec_x <- (z_2 - z_1) / 2
  vec_y <- (z_1 - z_2) / 2
  intersec_x <- z_1 + vec_x
  intersec_y <- z_2 + vec_y
  prob_loss <- sqrt((intersec_x^2 + intersec_y^2)) / sqrt(2)
  return(prob_loss)
}


#' Generates simulated data according to causal model with binary treatment
#' The standart of care is defined as treatment 1 if z_1 > 0 otherwise 0
#' @param seed Intager, seed.
#' @param n Intager, number of simulated individuals.
#' @param f decision function for counterfactual outcome
#' @examples
#' generate_data_binary(seed = 2, n = 100, f = descision_function)
generate_data <- function(n = 200, fun) {
  x <- runif(n, min = 0, max = 1)
  z_1 <- runif(n, min = 0, max = 1)
  z_2 <- runif(n, min = 0, max = 1)
  treat <- rbinom(n, prob = z_1, size = 1) # nolint: line_length_linter.
  prob <-ifelse(treat == TRUE, lossfunction_treat(z_1, z_2), 1 - lossfunction_treat(z_1, z_2)) # nolint: line_length_linter,    
  #prob <- ifelse(z_1 + z_2 > 1, (z_1^2 + z_2^2) / 4, expit(1 - z_1 + z_2 + treat * (z_1 + z_2))) # nolint: line_length_linter, commented_code_linter.
  y <- rbinom(n, prob = prob, size = 1)
  std_care_avg <- mean(prob)

  #simulate counterfactual data based on ground truth
  treat_c <- fun(data.frame(z_1, z_2))
  prob_c <-ifelse(treat_c == TRUE, lossfunction_treat(z_1, z_2), 1 - lossfunction_treat(z_1, z_2)) # nolint: line_length_linter,    
  #prob_c <- ifelse(z_1 + z_2 > 1, (z_1^2 + z_2 ^ 2) / 4, expit(1 - z_1 + z_2 + treat_c * (z_1 + z_2))) # nolint: line_length_linter, commented_code_linter.
  y_c <- rbinom(n, prob = prob_c, size = 1)
  countfact_avg <- mean(prob_c)

  sim_data <- data.frame(x, z_1, z_2, treat, y, treat_c, y_c) # nolint: line_length_linter.
  return(list(sim_data, std_care_avg, countfact_avg))
}


#' calculate naive IPW expected value
#'
#' @param data A list of vectors that cointains Z,T,Y,X.
#' @param f A decision function that take Z as an imput.
#' @return Expected value of Y with interfention T=f(Z) .
#' @examples
#' calculate_naive_ipw(df, decision_function)
calculate_naive_ipw <- function(df, f, plot = FALSE, estimated_weights = TRUE) {
  subindex <- f(df) == df$treat

  if (estimated_weights) {
    model <-
      glm(subindex ~ z_1 + z_2 + z_1 * z_2, family = binomial(link = "logit"), data = df) # nolint
    weights <-
      predict(model, df[subindex, c("z_1", "z_2")], type = "response")
  }else {
    weights <- ifelse(df$z_1 + df$z_2 < 1, df$z_1, 1 - df$z_1)
    weights <- weights[subindex]
  }
  expected_y <- sum(df$y[subindex] / weights) / nrow(df)

  if (plot) {
    plot(z_1 ~ z_2,
         xlab = "z2",
         ylab = "z1",
         main = "weights of used datapoints",
         col = ifelse(df$y, "#FFC20A", "#0C7BDC"),
         pch = ifelse(df$treat, 3, 1),
         cex = 1 / weights,
         data = df[subindex, ])
  }
  return(expected_y)
}

#' calculate expected value based on IPW
#'
#' @param data A list of vectors that cointains z,treat,y,x
#' @param f A decision function that take Z as an imput.
#' @return Expected value of Y with interfention treat=f(z) .
#' @examples
#' calculate_naive_ipw(df, decision_function)
calculate_ipw <- function(df, f, plot = FALSE, estimated_weights = TRUE) {
  subindex <- f(df) == df$treat
  if (estimated_weights) {
    model <-
      glm(df$treat ~ z_1 + z_2 + z_1 * z_2, family = binomial(link = "logit"), data = df) # nolint
    p_hat <-
      predict(model, df[subindex, c("z_1", "z_2")], type = "response")
    weights <- df$treat[subindex] / p_hat + (1 - df$treat[subindex]) / (1 - p_hat) # nolint: line_length_linter.

  }else {
    weights <- df$treat[subindex] / df$z_1[subindex] + (1 - df$treat[subindex]) / (1 - df$z_1[subindex]) # nolint: line_length_linter.
  }


  expected_y <- sum(df$y[subindex] * weights) / nrow(df)

  if (plot) {
    plot(z_1 ~ z_2,
         xlab = "z2",
         ylab = "z1",
         main = "weights of used datapoints",
         col = ifelse(df$y, "#FFC20A", "#0C7BDC"),
         pch = ifelse(df$treat, 3, 1),
         cex = weights,
         data = df[subindex, ])
  }

  return(expected_y)
}

#various types of g-computations
#of the form P(y=1 | f(z), z )
calculate_naive_g_comp <- function(df, f) {
  model <-
    glm(df$y ~ z_1 + z_2 + z_1 * z_2 + treat  + z_1 * treat + z_2 * treat + z_1 * z_2 * treat,  # nolint: line_length_linter.
        family = binomial(link = "logit"),
        data = df)

  df$treat <- as.integer(f(df))

  weights <-
    predict(model, df[c("z_1", "z_2", "treat")], type = "response")

  expected_y <- sum(weights) / nrow(df)
  return(expected_y)
}

calculate_g_comp <- function(df, f, estimated_weights = TRUE) {
  subindex <- f(df) == df$treat

  if (estimated_weights) {
    model <-
      glm(df$y ~ z_1 + z_2 + z_1 * z_2 + treat + z_1 * treat + z_2 * treat + z_1 * z_2 * treat,  # nolint: line_length_linter.
          family = binomial(link = "logit"),
          data = df)

    p_hats <-
      predict(model, df[subindex, c("z_1", "z_2", "treat")], type = "response") # nolint: line_length_linter.
    weights <-  df$y[subindex] * p_hats + (1 - df$y[subindex]) * (1 - p_hats)

    hist(weights)
  } else {
    p_hats <-ifelse(df$treat[subindex]== TRUE, lossfunction_treat(df$z_1, df$z_2), 1 - lossfunction_treat(df$z_1, df$z_2)) # nolint: line_length_linter,    
    weights <-  df$y[subindex] * p_hats + (1 - df$y[subindex]) * (1 - p_hats)
    hist(weights)
  }
  expected_y <- sum(weights) / nrow(df)
  return(expected_y)
}

#' calculate expected value based double robust estimator
#'
#' @param data A list of vectors that cointains z,treat,y,x
#' @param f A decision function that take Z as an imput.
#' @return Expected value of Y with interfention treat=f(z) .
#' @examples
#' calculate_naive_ipw(df, decision_function)
calculate_dr <- function(df, f, plot = FALSE, estimated_weights = TRUE) {
  subindex <- f(df) == df$treat
  if (estimated_weights) {
    model_ipw <-
      glm(df$treat ~ z_1 + z_2 + z_1 * z_2, family = binomial(link = "logit"), data = df) # nolint
    p_hat <-
      predict(model_ipw, df[subindex, c("z_1", "z_2")], type = "response")
    weights_ipw <- df$treat[subindex] / p_hat + (1 - df$treat[subindex]) / (1 - p_hat) # nolint: line_length_linter.

  }else {
    weights_ipw <- df$treat[subindex] / df$z_1[subindex] + (1 - df$treat[subindex]) / (1 - df$z_1[subindex]) # nolint: line_length_linter.
  }

  model_m <-
    glm(df$y ~ z_1 + z_2 + z_1 * z_2 + treat + z_1 * treat + z_2 * treat + z_1 * z_2 * treat,  # nolint: line_length_linter.
        family = binomial(link = "logit"),
        data = df)

  df$treat <- as.integer(f(df))

  y_hat <-
    predict(model_m, df[c("z_1", "z_2", "treat")], type = "response")

  expected_y <- (sum((df$y[subindex] - y_hat[subindex]) * weights_ipw) + sum(y_hat)) / nrow(df) # nolint: line_length_linter.

  return(expected_y)
}