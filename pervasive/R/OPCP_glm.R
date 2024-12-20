#' Calculate Observed Proportion of Concordant Pairs (OPCP)
#'
#' This function calculates the Observed Proportion of Concordant Pairs (OPCP) using Kendall's Tau as a measure of association. The pervasive functions also provide the OPCP.
#'
#' @param formula A formula specifying the dependent and independent variables.
#' @param data A data frame containing the variables specified in the formula.
#'
#' @return A numeric value representing the OPCP.
#' @examples
#' Example using the spi dataset from the psychTools package
#' library(psychTools)
#' sc <- psych::scoreVeryFast(spi.keys, spi)
#' spi_sc <- cbind(spi, sc)
#' spi_sc_vars <- spi_sc |>
#'   dplyr::select(sex, Agree, Consc, Neuro, Extra, Open)
#'   spi_sc_vars$sex = spi_sc_vars$sex -1
#'
#' formula <- sex ~ Agree + Consc + Neuro + Extra + Open
#' OPCP_glm(formula = formula, data = spi_sc_vars)
#'
#' @export
OPCP_glm <- function(formula, data) {
  # Step 1: Parse formula to extract dependent and independent variables
  dependent_var <- all.vars(formula)[1]
  independent_vars <- all.vars(formula)[-1]

  data <- na.omit(data[, c(dependent_var, independent_vars)])

  # Step 2: Run regression model and extract predictions
  model <- glm(formula, family=binomial(link='logit'), data = data)
  predictions <- predict(model, data)

  # Step 3: Calculate Kendall's Tau and OPCP
  observed <- data[[dependent_var]]
  tau <- cor(predictions, observed, method = "kendall", use = "pairwise.complete.obs")
  OPCP <- tau / 2 + 0.5

  return(OPCP)
}
