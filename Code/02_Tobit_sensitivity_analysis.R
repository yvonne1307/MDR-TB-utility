# ================================================================
# MDR-TB health utility analysis
# 1. Primary right-censored Tobit model
# 2. OLS with HC3 robust standard errors (same covariates as primary)
# 3. Treatment-status sensitivity Tobit model
# 4. Age-group sensitivity Tobit model
# ================================================================

library(readxl)
library(car)
library(censReg)
library(lmtest)
library(sandwich)
library(miscTools)

# ----------------------------
# 0. Data preparation
# ----------------------------

data_path <- paste0(
  "C:/Users/patients characteristics 0923.xlsx"
)

output_dir <- file.path(getwd(), "Tobit_analysis_outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read_excel(data_path)


# Continuous age, interpreted per 10-year increase.
dat$age10 <- dat$age / 10

# Convert categorical variables to unordered factors.
categorical_variables <- c(
  "sex",
  "occupation",
  "education_level",
  "annual_household_income",
  "medical_insurance_type",
  "self-reported_treatment_status",
  "self-reported_treatment_response",
  "adverse_reaction",
  "perceived_discrimination",
  "original_age_groups"
)

dat[categorical_variables] <- lapply(
  dat[categorical_variables],
  factor,
  ordered = FALSE
)

# Reference categories.
dat$sex <- relevel(dat$sex, ref = "Male")
dat$occupation <- relevel(
  dat$occupation,
  ref = "Farmer/Worker/Service worker"
)
dat$education_level <- relevel(
  dat$education_level,
  ref = "middle School or below"
)
dat$annual_household_income <- relevel(
  dat$annual_household_income,
  ref = "＜30,000 RMB"
)
dat$medical_insurance_type <- relevel(
  dat$medical_insurance_type,
  ref = paste0("The basic medical insurance system for urban and rural residents")
)
dat$`self-reported_treatment_status` <- relevel(
  dat$`self-reported_treatment_status`,
  ref = "Ongoing")

dat$`self-reported_treatment_response` <- relevel(
  dat$`self-reported_treatment_response`,
  ref = "Improved"
)
dat$adverse_reaction <- relevel(
  dat$adverse_reaction,
  ref = "No"
)
dat$perceived_discrimination <- relevel(
  dat$perceived_discrimination,
  ref = "No"
)
dat$original_age_groups <- relevel(
  dat$original_age_groups,
  ref = "18-25"
)

# ----------------------------
# 1. Reusable result functions
# ----------------------------

# Extract Tobit beta estimates, conventional maximum-likelihood standard
# errors, Wald 95% confidence intervals, and two-sided exact p values.
extract_tobit_coefficients <- function(model, analysis_name) {
  model_beta <- coef(model)
  model_vcov <- vcov(model)
  model_se <- sqrt(diag(model_vcov))
  model_z <- model_beta / model_se
  model_p <- 2 * pnorm(abs(model_z), lower.tail = FALSE)
  critical_value <- qnorm(0.975)

  result <- data.frame(
    Analysis = analysis_name,
    Variable = names(model_beta),
    Estimate = as.numeric(model_beta),
    Standard_error = as.numeric(model_se),
    CI_95_lower = as.numeric(model_beta - critical_value * model_se),
    CI_95_upper = as.numeric(model_beta + critical_value * model_se),
    Z_value = as.numeric(model_z),
    P_value = as.numeric(model_p),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # logSigma is a distributional parameter rather than a covariate effect.
  # Keep it in the complete output but flag it so that it can be omitted from
  # the publication regression table.
  result$Parameter_type <- ifelse(
    grepl("logSigma|sigma", result$Variable, ignore.case = TRUE),
    "Distribution parameter",
    "Regression coefficient"
  )

  result$Estimate_display <- sprintf("%.3f", result$Estimate)
  result$SE_display <- sprintf("%.3f", result$Standard_error)
  result$CI_display <- sprintf(
    "%.3f to %.3f",
    result$CI_95_lower,
    result$CI_95_upper
  )
  result$P_display <- ifelse(
    result$P_value < 0.001,
    "<0.001",
    sprintf("%.3f", result$P_value)
  )

  result
}

# Extract log-likelihood, AIC, parameter count, and censoring counts for a
# right-censored Tobit model. Counts are based on the complete cases actually
# available for the corresponding model formula.
extract_tobit_fit <- function(
  model,
  model_formula,
  analysis_data,
  analysis_name,
  upper_limit = 1
) {
  model_frame <- model.frame(
    model_formula,
    data = analysis_data,
    na.action = na.omit
  )
  outcome_used <- model.response(model_frame)
  model_loglik <- logLik(model)

  data.frame(
    Analysis = analysis_name,
    N = length(outcome_used),
    Uncensored_N = sum(outcome_used < upper_limit),
    Right_censored_N = sum(outcome_used >= upper_limit),
    Right_censored_percent = 100 *
      sum(outcome_used >= upper_limit) / length(outcome_used),
    Log_likelihood = as.numeric(model_loglik),
    Number_of_parameters = attr(model_loglik, "df"),
    AIC = AIC(model),
    stringsAsFactors = FALSE
  )
}

# Standardize the output from car::vif(). When a model contains a
# multi-category factor, car::vif() returns GVIF, Df, and adjusted GVIF.
# If every term has one degree of freedom, it returns ordinary VIF values.
extract_vif_gvif <- function(lm_model, analysis_name) {
  vif_object <- car::vif(lm_model)

  if (is.matrix(vif_object)) {
    result <- data.frame(
      Analysis = analysis_name,
      Variable = rownames(vif_object),
      GVIF = as.numeric(vif_object[, "GVIF"]),
      Df = as.numeric(vif_object[, "Df"]),
      GVIF_adjusted = as.numeric(
        vif_object[, "GVIF^(1/(2*Df))"]
      ),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  } else {
    result <- data.frame(
      Analysis = analysis_name,
      Variable = names(vif_object),
      GVIF = as.numeric(vif_object),
      Df = 1,
      GVIF_adjusted = sqrt(as.numeric(vif_object)),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }

  result
}

# Extract OLS estimates using HC3 heteroskedasticity-robust standard errors.
# Confidence intervals and p values use the residual t distribution.
extract_ols_hc3_coefficients <- function(
  ols_model,
  analysis_name
) {
  hc3_vcov <- sandwich::vcovHC(ols_model, type = "HC3")
  hc3_test <- lmtest::coeftest(
    ols_model,
    vcov. = hc3_vcov
  )
  residual_df <- df.residual(ols_model)
  critical_value <- qt(0.975, df = residual_df)

  result <- data.frame(
    Analysis = analysis_name,
    Variable = rownames(hc3_test),
    Estimate = as.numeric(hc3_test[, 1]),
    Standard_error_HC3 = as.numeric(hc3_test[, 2]),
    CI_95_lower = as.numeric(
      hc3_test[, 1] - critical_value * hc3_test[, 2]
    ),
    CI_95_upper = as.numeric(
      hc3_test[, 1] + critical_value * hc3_test[, 2]
    ),
    T_value = as.numeric(hc3_test[, 3]),
    P_value = as.numeric(hc3_test[, 4]),
    Residual_df = residual_df,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  result$Estimate_display <- sprintf("%.3f", result$Estimate)
  result$SE_display <- sprintf("%.3f", result$Standard_error_HC3)
  result$CI_display <- sprintf(
    "%.3f to %.3f",
    result$CI_95_lower,
    result$CI_95_upper
  )
  result$P_display <- ifelse(
    result$P_value < 0.001,
    "<0.001",
    sprintf("%.3f", result$P_value)
  )

  result
}

# Extract marginal effects of the primary Tobit model on the expected
# observed utility value. censReg provides the S3 method, whereas the margEff
# generic is exported by miscTools. By default these effects are evaluated at
# the covariate means; they are therefore marginal effects at the means rather
# than sample-average marginal effects.
extract_tobit_marginal_effects <- function(
  model,
  analysis_name
) {
  marginal_effect_object <- miscTools::margEff(
    model,
    calcVCov = TRUE
  )

  marginal_effect_matrix <- summary(
    marginal_effect_object
  )

  marginal_effect_df <- attr(
    marginal_effect_object,
    "df.residual"
  )

  if (
    is.null(marginal_effect_df) ||
      length(marginal_effect_df) == 0 ||
      !is.finite(marginal_effect_df)
  ) {
    critical_value <- qnorm(0.975)
  } else {
    critical_value <- qt(
      0.975,
      df = marginal_effect_df
    )
  }

  result <- data.frame(
    Analysis = analysis_name,
    Variable = rownames(marginal_effect_matrix),
    Adjusted_marginal_effect = as.numeric(
      marginal_effect_matrix[, 1]
    ),
    Standard_error = as.numeric(
      marginal_effect_matrix[, 2]
    ),
    Test_statistic = as.numeric(
      marginal_effect_matrix[, 3]
    ),
    P_value = as.numeric(
      marginal_effect_matrix[, 4]
    ),
    Evaluation = "At covariate means",
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  result$CI_95_lower <-
    result$Adjusted_marginal_effect -
    critical_value * result$Standard_error

  result$CI_95_upper <-
    result$Adjusted_marginal_effect +
    critical_value * result$Standard_error

  # The intercept does not represent a covariate marginal effect.
  result <- result[
    result$Variable != "(Intercept)",
  ]

  result$Marginal_effect_display <- sprintf(
    "%.3f",
    result$Adjusted_marginal_effect
  )
  result$SE_display <- sprintf(
    "%.3f",
    result$Standard_error
  )
  result$CI_display <- sprintf(
    "%.3f to %.3f",
    result$CI_95_lower,
    result$CI_95_upper
  )
  result$P_display <- ifelse(
    result$P_value < 0.001,
    "<0.001",
    sprintf("%.3f", result$P_value)
  )

  result
}

# ----------------------------
# 2. Targeted treatment-variable collinearity assessment
# ----------------------------

treatment_crosstab <- table(
  Treatment_status = dat$`self-reported_treatment_status`,
  Treatment_response = dat$`self-reported_treatment_response`
)

treatment_row_percent <- prop.table(
  treatment_crosstab,
  margin = 1
) * 100

treatment_chisq <- chisq.test(treatment_crosstab)

treatment_cramers_v <- sqrt(
  as.numeric(treatment_chisq$statistic) /
    (
      sum(treatment_crosstab) *
        min(
          nrow(treatment_crosstab) - 1,
          ncol(treatment_crosstab) - 1
        )
    )
)

treatment_association_results <- data.frame(
  Pearson_chi_square = as.numeric(treatment_chisq$statistic),
  Degrees_of_freedom = as.numeric(treatment_chisq$parameter),
  P_value = treatment_chisq$p.value,
  Cramers_V = treatment_cramers_v,
  Any_expected_count_below_5 = any(treatment_chisq$expected < 5)
)

# Optional candidate model containing both treatment variables. This is used
# only to demonstrate their overlap and is not a final Tobit model.
candidate_collinearity_formula <- Utility ~
  age10 +
  sex +
  occupation +
  education_level +
  annual_household_income +
  medical_insurance_type +
  `self-reported_treatment_status` +
  `self-reported_treatment_response` +
  adverse_reaction +
  perceived_discrimination

candidate_collinearity_lm <- lm(
  candidate_collinearity_formula,
  data = dat
)

candidate_vif_gvif_table <- extract_vif_gvif(
  candidate_collinearity_lm,
  "Candidate model containing both treatment variables"
)

# ----------------------------
# 3. Primary Tobit analysis
# ----------------------------

primary_tobit_formula <- Utility ~
  age10 +
  sex +
  occupation +
  education_level +
  annual_household_income +
  medical_insurance_type +
  `self-reported_treatment_response` +
  adverse_reaction +
  perceived_discrimination

# VIF/GVIF for the final primary model.
primary_vif_lm <- lm(
  primary_tobit_formula,
  data = dat
)

primary_vif_gvif_table <- extract_vif_gvif(
  primary_vif_lm,
  "Primary Tobit model"
)

# Right-censored Tobit model: values equal to 1 are right-censored; all
# values below 1, including negative utility values, are uncensored.
primary_tobit_model <- censReg(
  primary_tobit_formula,
  left = -Inf,
  right = 1,
  data = dat
)

primary_tobit_coefficient_table <- extract_tobit_coefficients(
  primary_tobit_model,
  "Primary Tobit model"
)

primary_tobit_fit_table <- extract_tobit_fit(
  primary_tobit_model,
  primary_tobit_formula,
  dat,
  "Primary Tobit model"
)

# Publication table excludes the distribution parameter logSigma.
primary_tobit_publication_table <- subset(
  primary_tobit_coefficient_table,
  Parameter_type == "Regression coefficient"
)

# Adjusted marginal effects on the expected observed utility value. These are
# evaluated at the covariate means and aid interpretation in absolute utility
# units. They should not be labelled as average marginal effects.
primary_tobit_marginal_effect_table <-
  extract_tobit_marginal_effects(
    primary_tobit_model,
    "Primary Tobit model"
  )

# ----------------------------
# 4. OLS + HC3 robustness analysis
# ----------------------------

# This model uses exactly the same outcome and covariates as the primary
# Tobit model.
primary_ols_hc3_model <- lm(
  primary_tobit_formula,
  data = dat
)

primary_ols_hc3_coefficient_table <-
  extract_ols_hc3_coefficients(
    primary_ols_hc3_model,
    "OLS with HC3 robust standard errors"
  )

# ----------------------------
# 5. Treatment-status sensitivity Tobit analysis
# ----------------------------

# Treatment status replaces self-reported treatment response. All other
# covariates are kept identical to the primary model.
treatment_status_sensitivity_formula <- Utility ~
  age10 +
  sex +
  occupation +
  education_level +
  annual_household_income +
  medical_insurance_type +
  `self-reported_treatment_status` +
  adverse_reaction +
  perceived_discrimination

treatment_status_sensitivity_vif_lm <- lm(
  treatment_status_sensitivity_formula,
  data = dat
)

treatment_status_sensitivity_vif_gvif_table <-
  extract_vif_gvif(
    treatment_status_sensitivity_vif_lm,
    "Treatment-status sensitivity Tobit model"
  )

treatment_status_sensitivity_tobit_model <- censReg(
  treatment_status_sensitivity_formula,
  left = -Inf,
  right = 1,
  data = dat
)

treatment_status_sensitivity_coefficient_table <-
  extract_tobit_coefficients(
    treatment_status_sensitivity_tobit_model,
    "Treatment-status sensitivity Tobit model"
  )

treatment_status_sensitivity_fit_table <- extract_tobit_fit(
  treatment_status_sensitivity_tobit_model,
  treatment_status_sensitivity_formula,
  dat,
  "Treatment-status sensitivity Tobit model"
)

treatment_status_sensitivity_publication_table <- subset(
  treatment_status_sensitivity_coefficient_table,
  Parameter_type == "Regression coefficient"
)

# ----------------------------
# 6. Age-group sensitivity Tobit analysis
# ----------------------------

# Original age groups replace continuous age. All other covariates are kept
# identical to the primary model.
age_group_sensitivity_formula <- Utility ~
  original_age_groups +
  sex +
  occupation +
  education_level +
  annual_household_income +
  medical_insurance_type +
  `self-reported_treatment_response` +
  adverse_reaction +
  perceived_discrimination

age_group_sensitivity_vif_lm <- lm(
  age_group_sensitivity_formula,
  data = dat
)

age_group_sensitivity_vif_gvif_table <- extract_vif_gvif(
  age_group_sensitivity_vif_lm,
  "Age-group sensitivity Tobit model"
)

age_group_sensitivity_tobit_model <- censReg(
  age_group_sensitivity_formula,
  left = -Inf,
  right = 1,
  data = dat
)

age_group_sensitivity_coefficient_table <-
  extract_tobit_coefficients(
    age_group_sensitivity_tobit_model,
    "Age-group sensitivity Tobit model"
  )

age_group_sensitivity_fit_table <- extract_tobit_fit(
  age_group_sensitivity_tobit_model,
  age_group_sensitivity_formula,
  dat,
  "Age-group sensitivity Tobit model"
)

age_group_sensitivity_publication_table <- subset(
  age_group_sensitivity_coefficient_table,
  Parameter_type == "Regression coefficient"
)

# ----------------------------
# 7. Combined review tables
# ----------------------------

all_tobit_fit_statistics <- rbind(
  primary_tobit_fit_table,
  treatment_status_sensitivity_fit_table,
  age_group_sensitivity_fit_table
)

all_final_model_vif_gvif <- rbind(
  primary_vif_gvif_table,
  treatment_status_sensitivity_vif_gvif_table,
  age_group_sensitivity_vif_gvif_table
)

# Compact coefficient tables intended for manuscript/supplement preparation.
select_publication_columns <- function(x) {
  x[, c(
    "Analysis",
    "Variable",
    "Estimate",
    "Standard_error",
    "CI_95_lower",
    "CI_95_upper",
    "P_value",
    "Estimate_display",
    "SE_display",
    "CI_display",
    "P_display"
  )]
}

primary_tobit_publication_table <- select_publication_columns(
  primary_tobit_publication_table
)
treatment_status_sensitivity_publication_table <-
  select_publication_columns(
    treatment_status_sensitivity_publication_table
  )
age_group_sensitivity_publication_table <-
  select_publication_columns(
    age_group_sensitivity_publication_table
  )

# ----------------------------
# 8. Display results in the R console
# ----------------------------

primary_tobit_publication_table
primary_tobit_fit_table
primary_vif_gvif_table
primary_tobit_marginal_effect_table

primary_ols_hc3_coefficient_table

treatment_status_sensitivity_publication_table
treatment_status_sensitivity_fit_table
treatment_status_sensitivity_vif_gvif_table

age_group_sensitivity_publication_table
age_group_sensitivity_fit_table
age_group_sensitivity_vif_gvif_table

treatment_crosstab
treatment_row_percent
treatment_association_results

# ----------------------------
# 9. Export all requested results as CSV files
# ----------------------------

write.csv(
  primary_tobit_publication_table,
  file.path(output_dir, "01_primary_Tobit_coefficients.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  primary_tobit_fit_table,
  file.path(output_dir, "02_primary_Tobit_model_fit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  primary_vif_gvif_table,
  file.path(output_dir, "03_primary_Tobit_VIF_GVIF.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  primary_tobit_marginal_effect_table,
  file.path(
    output_dir,
    "04_primary_Tobit_adjusted_marginal_effects.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  primary_ols_hc3_coefficient_table,
  file.path(output_dir, "05_primary_OLS_HC3_coefficients.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  treatment_status_sensitivity_publication_table,
  file.path(
    output_dir,
    "06_treatment_status_sensitivity_Tobit_coefficients.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  treatment_status_sensitivity_fit_table,
  file.path(
    output_dir,
    "07_treatment_status_sensitivity_Tobit_model_fit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  treatment_status_sensitivity_vif_gvif_table,
  file.path(
    output_dir,
    "08_treatment_status_sensitivity_VIF_GVIF.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  age_group_sensitivity_publication_table,
  file.path(
    output_dir,
    "09_age_group_sensitivity_Tobit_coefficients.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  age_group_sensitivity_fit_table,
  file.path(
    output_dir,
    "10_age_group_sensitivity_Tobit_model_fit.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  age_group_sensitivity_vif_gvif_table,
  file.path(
    output_dir,
    "11_age_group_sensitivity_VIF_GVIF.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  all_tobit_fit_statistics,
  file.path(output_dir, "12_all_Tobit_model_fit_statistics.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  all_final_model_vif_gvif,
  file.path(output_dir, "13_all_final_models_VIF_GVIF.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  candidate_vif_gvif_table,
  file.path(
    output_dir,
    "14_candidate_model_with_both_treatment_variables_VIF_GVIF.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  as.data.frame.matrix(treatment_crosstab),
  file.path(output_dir, "15_treatment_status_response_crosstab.csv"),
  row.names = TRUE,
  fileEncoding = "UTF-8"
)
write.csv(
  treatment_association_results,
  file.path(output_dir, "16_treatment_variable_association.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message(
  "Analysis complete. Results were saved to: ",
  normalizePath(output_dir, winslash = "/", mustWork = FALSE)
)
