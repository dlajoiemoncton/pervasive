#' Calculate Observed Proportion of Concordant Pairs (OPCP)
#'
#' This function provides a matrix (like a correlation matrix) of Observed Proportions of Concordant Pairs (OPCPs) .
#'
#' @param data A data frame containing the variables specified in the formula.
#'
#' @return A matrix of OPCPs.
#' @examples
#' # Example using the spi dataset from the psychTools package
#' library(psychTools)
#' sc <- psych::scoreVeryFast(spi.keys, spi)
#' spi_sc <- cbind(spi, sc)
#' spi_sc_age_sex_B5 <- spi_sc |>
#'   dplyr::select(age, sex, Agree, Consc, Neuro, Extra, Open)
#'
#'
#' OPCP_mat(data = spi_sc_age_sex_B5)
#'
#' @export
OPCP_mat <- function(data) {

  # Compute the Kendall's tau correlation matrix
  tau_matrix <- cor(data, method = "kendall", use = "pairwise.complete.obs")

  # Apply the OPCP transformation
  opcp_mat <- (tau_matrix / 2) + 0.5

  # Return the OPCP matrix
  return(opcp_mat)

}
