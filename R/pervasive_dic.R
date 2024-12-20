#' Association Rule Mining With Dichotomized Data
#'
#' This function provides the same information as pervasive_tric but with dichotomized data instead of trichotomized data. Trichotomized data is generally to be preferred when possible. Researchers interested in mixing trichotomized variables with dichotomized variables should use the arules package directly.
#'
#' @param formula A formula specifying the dependent and independent variables.
#' @param data A data frame containing the variables specified in the formula.
#'
#' @return @return A list with the following components:
#' \itemize{
#'   \item \code{OPCP}: Observed proportion of concordant pairs.
#'   \item \code{adj_r_squared}: Adjusted R-squared value for the regression model.
#'   \item \code{exact_match_lhs}, \code{exact_match_rhs}: The left and right-hand side of the rule suggested by the regression model, respectively
#'   \item \code{exact_match_quality}: Quality metrics for the rule suggested by the regression.
#'   \item \code{exact_match_lhs_opp}, \code{exact_match_rhs_opp}: The left and right-hand side of the rule suggested by the low end of the regression model, respectively
#'   \item \code{exact_match_quality_opp}: Quality metrics for the rule suggested by the low end of the regression.
#'   \item \code{top_rule_lhs}, \code{top_rule_rhs}, \code{top_rule_quality}: Information relevant to the highest lift rule meeting min_support for high values of the dependent variable.
#'   \item \code{top_rule_opp_lhs}, \code{top_rule_opp_rhs}, \code{top_rule_opp_quality}: Information relevant to the highest lift rule meeting min_support for low values of the dependent variable.
#'   \item \code{quality_table}: A table summarizing the quality statistics for extracted association rules.
#'   \item \code{freq_tables}: Frequency tables (cutoffs and membership) for dichotomization binning.
#' }
#' @examples
#' # Example using the spi dataset from the psychTools package
#' library(psychTools)
#' sc <- psych::scoreVeryFast(spi.keys, spi)
#' spi_sc <- cbind(spi, sc)
#' spi_sc_vars <- spi_sc |>
#'   dplyr::select(age, Agree, Consc, Neuro, Extra, Open)
#'
#' formula <- age ~ Agree + Consc + Neuro + Extra + Open
#' example <- pervasive_dic(formula = formula, data = spi_sc_vars)
#' example
#'
#' @import arules
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tibble rownames_to_column
#'
#'
#'
#' @export
pervasive_dic <- function(formula, data, min_support) {
  # Step 1: Parse formula to extract dependent and independent variables
  dependent_var <- all.vars(formula)[1]  # Dependent variable
  independent_vars <- all.vars(formula)[-1]  # Independent variables

  data <- na.omit(data[, c(dependent_var, independent_vars)])

  # Step 2: Run regression model and extract predictions
  model <- lm(formula, data = data)
  predictions <- predict(model, data)

  # Step 3: Calculate Kendall's Tau and OPCP
  observed <- data[[dependent_var]]
  tau <- cor(predictions, observed, method = "kendall", use = "pairwise.complete.obs")
  OPCP <- tau / 2 + 0.5

  # Step 4: Extract adjusted R-squared
  adj_r_squared <- summary(model)$adj.r.squared

  # Step 5: Identify significant variables and their valence
  significant_vars <- summary(model)$coefficients |>
    as.data.frame() |>
    tibble::rownames_to_column("Variable") |>
    dplyr::filter(`Pr(>|t|)` < 0.05) |>
    dplyr::filter(Variable != "(Intercept)") |>
    dplyr::mutate(Valence = ifelse(Estimate > 0, "high", "low")) |>
    dplyr::mutate(OppValence = ifelse(Estimate > 0, "low", "high"))

  # Step 6: Dichotomize dependent and independent variables
  dichotomize <- function(x) {
    if (length(unique(x)) <= 2) {
      # Already binary; relabel as "low" and "high"
      factor(x, levels = c(0, 1), labels = c("low", "high"), ordered = TRUE)
    } else {
    discretize(x, breaks = 2, labels = c("low", "high"),
               include.lowest = TRUE, ordered_result = TRUE)
  }}

  dichotomized_dep <- data |>
    dplyr::mutate(across(all_of(dependent_var), dichotomize, .names = "dich_{col}"))
  transactions_dep <- as(dichotomized_dep |> select(starts_with("dich_")), "transactions")

  dichotomized_ind <- data |>
    dplyr::mutate(across(all_of(independent_vars), dichotomize, .names = "dich_{col}"))
  transactions_ind <- as(dichotomized_ind |> select(starts_with("dich_")), "transactions")

  transactions <- as(cbind(dichotomized_dep, dichotomized_ind) |>
                       select(starts_with("dich_")), "transactions")

  # Step 7: Generate concatenated strings for variables
  significant_vars <- significant_vars |>
    dplyr::mutate(dich_string = paste0("dich_", Variable, "=", Valence)) |>
    dplyr::mutate(dich_string_opp = paste0("dich_", Variable, "=", OppValence))

  dep_var_high <- paste0("dich_", dependent_var, "=high")
  dep_var_low <- paste0("dich_", dependent_var, "=low")

  dich_vars <- significant_vars$dich_string
  dich_vars_opp <- significant_vars$dich_string_opp

  # Step 8: Extract relevant association rules using apriori
  rules <- apriori(transactions, support = 1 / nrow(transactions), confidence = 0,
                   appearance = list(lhs = colnames(transactions_ind),
                                     rhs = colnames(transactions_dep)))

  rules_w_all_lhs <- subset(rules, lhs %ain% dich_vars)
  rules_w_only_lhs <- subset(rules_w_all_lhs, lhs %oin% dich_vars)
  exact_match_rule <- subset(rules_w_only_lhs, rhs %in% dep_var_high)

  rules_w_all_lhs_opp <- subset(rules, lhs %ain% dich_vars_opp)
  rules_w_only_lhs_opp <- subset(rules_w_all_lhs_opp, lhs %oin% dich_vars_opp)
  exact_match_rule_opp <- subset(rules_w_only_lhs_opp, rhs %in% dep_var_low)

  rules_with_dep_var_high <- subset(rules,
                                    rhs %pin% dep_var_high &
                                      quality(rules)$support >= min_support)
  rules_with_dep_var_low <- subset(rules,
                                   rhs %pin% dep_var_low &
                                     quality(rules)$support >= min_support)

  # Step 9: Extract specific rules and their qualities

  exact_match_lhs <- labels(lhs(exact_match_rule))
  exact_match_rhs <- labels(rhs(exact_match_rule))
  exact_match_quality <- quality(exact_match_rule)

  exact_match_lhs_opp <- labels(lhs(exact_match_rule_opp))
  exact_match_rhs_opp <- labels(rhs(exact_match_rule_opp))
  exact_match_quality_opp <- quality(exact_match_rule_opp)

  top_rule <- head(rules_with_dep_var_high, n = 1, by = "lift")
  top_rule_lhs <- labels(lhs(top_rule))
  top_rule_rhs <- labels(rhs(top_rule))
  top_rule_quality <- quality(top_rule)

  top_rule_opp <- head(rules_with_dep_var_low, n = 1, by = "lift")
  top_rule_opp_lhs <- labels(lhs(top_rule_opp))
  top_rule_opp_rhs <- labels(rhs(top_rule_opp))
  top_rule_opp_quality <- quality(top_rule_opp)


  # Helper function to safely extract rule quality metrics
  safe_extract <- function(quality_obj, metric) {
    if (!is.null(quality_obj) && length(quality_obj[[metric]]) > 0) {
      quality_obj[[metric]]
    } else {
      NA
    }
  }

  # Step 10: Create summary table for the rules
  quality_table <- data.frame(
    Description = c(
      "The rule suggested by your regression",
      "The rule suggested by the low end of your regression",
      "The rule with the highest lift for high values of your outcome",
      "The rule with the highest lift for low values of your outcome"
    ),
    LHS = c(
      if (!is.null(exact_match_rule) && length(exact_match_rule) > 0) {
        paste(labels(lhs(exact_match_rule)), collapse = " & ")
      } else { NA },
      if (!is.null(exact_match_rule_opp) && length(exact_match_rule_opp) > 0) {
        paste(labels(lhs(exact_match_rule_opp)), collapse = " & ")
      } else { NA },
      if (!is.null(top_rule) && length(top_rule) > 0) {
        paste(labels(lhs(top_rule)), collapse = " & ")
      } else { NA },
      if (!is.null(top_rule_opp) && length(top_rule_opp) > 0) {
        paste(labels(lhs(top_rule_opp)), collapse = " & ")
      } else { NA }
    ),
    RHS = c(
      if (!is.null(exact_match_rule) && length(exact_match_rule) > 0) {
        paste(labels(rhs(exact_match_rule)), collapse = " & ")
      } else { NA },
      if (!is.null(exact_match_rule_opp) && length(exact_match_rule_opp) > 0) {
        paste(labels(rhs(exact_match_rule_opp)), collapse = " & ")
      } else { NA },
      if (!is.null(top_rule) && length(top_rule) > 0) {
        paste(labels(rhs(top_rule)), collapse = " & ")
      } else { NA },
      if (!is.null(top_rule_opp) && length(top_rule_opp) > 0) {
        paste(labels(rhs(top_rule_opp)), collapse = " & ")
      } else { NA }
    ),
    Support = c(
      safe_extract(exact_match_quality, "support"),
      safe_extract(exact_match_quality_opp, "support"),
      safe_extract(quality(top_rule), "support"),
      safe_extract(quality(top_rule_opp), "support")
    ),
    Confidence = c(
      safe_extract(exact_match_quality, "confidence"),
      safe_extract(exact_match_quality_opp, "confidence"),
      safe_extract(quality(top_rule), "confidence"),
      safe_extract(quality(top_rule_opp), "confidence")
    ),
    Coverage = c(
      safe_extract(exact_match_quality, "coverage"),
      safe_extract(exact_match_quality_opp, "coverage"),
      safe_extract(top_rule_quality, "coverage"),
      safe_extract(top_rule_opp_quality, "coverage")
    ),
    Lift = c(
      safe_extract(exact_match_quality, "lift"),
      safe_extract(exact_match_quality_opp, "lift"),
      safe_extract(quality(top_rule), "lift"),
      safe_extract(quality(top_rule_opp), "lift")
    ),
    stringsAsFactors = FALSE
  )


  # Step 11: Generate frequency tables for dichotomized variables
  dic_freq <- function(x) {
    if (length(unique(x)) <= 2) {
      # Already binary;
      factor(x, levels = c(0, 1),  ordered = TRUE)
    } else {
    discretize(x, breaks = 2, include.lowest = TRUE, ordered_result = TRUE)
  }}

  dic_freq_binned <- data |>
    select(all_of(c(dependent_var, independent_vars))) |>
    dplyr::mutate(across(all_of(c(dependent_var, independent_vars)), dic_freq))

  freq_tables <- lapply(dic_freq_binned, table)


  # Step 12: Create final output
  result <- list(
    OPCP = OPCP,
    adj_r_squared = adj_r_squared,
    exact_match_lhs = exact_match_lhs,
    exact_match_rhs = exact_match_rhs,
    exact_match_quality = exact_match_quality,
    exact_match_lhs_opp = exact_match_lhs_opp,
    exact_match_rhs_opp = exact_match_rhs_opp,
    exact_match_quality_opp = exact_match_quality_opp,
    top_rule_lhs = top_rule_lhs,
    top_rule_rhs = top_rule_rhs,
    top_rule_quality = top_rule_quality,
    top_rule_opp_lhs = top_rule_opp_lhs,
    top_rule_opp_rhs = top_rule_opp_rhs,
    top_rule_opp_quality = top_rule_opp_quality,
    quality_table = quality_table,
    freq_tables = freq_tables
  )


  # Add class to result for custom printing
  class(result) <- "pervasive_dic"
  return(result)
}

#' @export
# Custom print function for pervasive object
print.pervasive_dic <- function(x, ...) {
  cat("\nThe adjusted R-square of your regression is",
      format(x$adj_r_squared, digits = 2),
      "with an observed percentage of concordant pairs (OPCP) of",
      format(x$OPCP * 100, digits = 4), "%.\n\n")

  print(x$quality_table, row.names = FALSE)
  invisible(x)
}
