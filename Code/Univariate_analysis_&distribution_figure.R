# ================================================================
# MDR-TB univariable analysis and utility distribution figure
# Data source: patients characteristics 0923.xlsx
# ================================================================

library(readxl)
library(ggplot2)

# ----------------------------
# 0. Paths and data
# ----------------------------

data_path <- paste0(
  "C:/Users/patients characteristics 0923.xlsx"
)

output_dir <- file.path(getwd(), "Univariate_analysis_outputs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read_excel(data_path)

required_columns <- c(
  "age",
  "original_age_groups",
  "sex",
  "education_level",
  "occupation",
  "annual_household_income",
  "medical_insurance_type",
  "self-reported_treatment_status",
  "self-reported_treatment_response",
  "adverse_reaction",
  "perceived_discrimination",
  "Utility"
)

missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0) {
  stop(
    "The following required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

dat$age <- as.numeric(dat$age)
dat$Utility <- as.numeric(dat$Utility)

# Set clinically meaningful display order for categorical variables.
dat$original_age_groups <- factor(
  dat$original_age_groups,
  levels = c("18-25", "26-35", "36-45", "46-55", "56-60")
)
dat$sex <- factor(
  dat$sex,
  levels = c("Male", "Female")
)
dat$education_level <- factor(
  dat$education_level,
  levels = c(
    "middle School or below",
    "Senior School or Technical Secondary School",
    "Junior college",
    "Bachelor's degree or above"
  )
)
dat$occupation <- factor(
  dat$occupation,
  levels = c(
    "Farmer/Worker/Service worker",
    "staff",
    "Self-employed",
    "Retiree/Unemployed",
    "Student",
    "Others"
  )
)
dat$annual_household_income <- factor(
  dat$annual_household_income,
  levels = c(
    "＜30,000 RMB",
    "30,000 - 80,000 RMB",
    "80,000 - 150,000 RMB",
    "＞150,000 RMB"
  )
)
dat$medical_insurance_type <- factor(
  dat$medical_insurance_type,
  levels = c(
    "The basic medical insurance system for urban and rural residents",
    "The basic medical insurance for urban employees",
    "Unclearly"
  )
)
dat$`self-reported_treatment_status` <- factor(
  dat$`self-reported_treatment_status`,
  levels = c("Ongoing", "Completed")
)
dat$`self-reported_treatment_response` <- factor(
  dat$`self-reported_treatment_response`,
  levels = c(
    "Improved",
    "Cured",
    "Unfavorable (no change/recurrence)"
  )
)
dat$adverse_reaction <- factor(
  dat$adverse_reaction,
  levels = c("No", "Yes")
)
dat$perceived_discrimination <- factor(
  dat$perceived_discrimination,
  levels = c("No", "Slight", "greater"),
  labels = c("No", "Slight", "Greater")
)

# ----------------------------
# 1. Overall utility summary
# ----------------------------

utility_complete <- dat$Utility[!is.na(dat$Utility)]

overall_utility_summary <- data.frame(
  N = length(utility_complete),
  Missing_N = sum(is.na(dat$Utility)),
  Mean = mean(utility_complete),
  Standard_deviation = sd(utility_complete),
  Median = median(utility_complete),
  Q1 = as.numeric(quantile(utility_complete, 0.25)),
  Q3 = as.numeric(quantile(utility_complete, 0.75)),
  Minimum = min(utility_complete),
  Maximum = max(utility_complete),
  Utility_equal_to_1_N = sum(utility_complete == 1),
  Utility_below_0_N = sum(utility_complete < 0),
  stringsAsFactors = FALSE
)

overall_utility_summary$Mean_SD_display <- sprintf(
  "%.3f ± %.3f",
  overall_utility_summary$Mean,
  overall_utility_summary$Standard_deviation
)
overall_utility_summary$Median_IQR_display <- sprintf(
  "%.3f (%.3f, %.3f)",
  overall_utility_summary$Median,
  overall_utility_summary$Q1,
  overall_utility_summary$Q3
)

# ----------------------------
# 2. Univariable categorical analyses
# ----------------------------

analysis_variables <- c(
  "original_age_groups",
  "sex",
  "education_level",
  "occupation",
  "annual_household_income",
  "medical_insurance_type",
  "self-reported_treatment_status",
  "self-reported_treatment_response",
  "adverse_reaction",
  "perceived_discrimination"
)

variable_labels <- c(
  original_age_groups = "Age group",
  sex = "Sex",
  education_level = "Education level",
  occupation = "Occupation",
  annual_household_income = "Annual household income",
  medical_insurance_type = "Medical insurance type",
  `self-reported_treatment_status` = "Self-reported treatment status",
  `self-reported_treatment_response` = "Self-reported treatment response",
  adverse_reaction = "Adverse reaction",
  perceived_discrimination = "Perceived discrimination"
)

analyse_one_categorical_variable <- function(
  data,
  variable_name,
  display_name
) {
  complete_index <- complete.cases(
    data[, c(variable_name, "Utility")]
  )
  analysis_data <- data[complete_index, c(variable_name, "Utility")]
  names(analysis_data) <- c("Group", "Utility")
  analysis_data$Group <- droplevels(factor(analysis_data$Group))

  group_levels <- levels(analysis_data$Group)
  number_groups <- length(group_levels)
  valid_n <- nrow(analysis_data)

  if (number_groups < 2) {
    stop(
      "Variable '", variable_name,
      "' has fewer than two observed categories."
    )
  }

  descriptive_rows <- do.call(
    rbind,
    lapply(group_levels, function(group_level) {
      group_utility <- analysis_data$Utility[
        analysis_data$Group == group_level
      ]

      data.frame(
        Variable = display_name,
        Variable_name = variable_name,
        Category = group_level,
        N = length(group_utility),
        Percent = 100 * length(group_utility) / valid_n,
        Mean = mean(group_utility),
        Standard_deviation = sd(group_utility),
        Median = median(group_utility),
        Q1 = as.numeric(quantile(group_utility, 0.25)),
        Q3 = as.numeric(quantile(group_utility, 0.75)),
        Valid_N = valid_n,
        Missing_N = nrow(data) - valid_n,
        stringsAsFactors = FALSE
      )
    })
  )

  descriptive_rows$Percent_display <- sprintf(
    "%.1f%%",
    descriptive_rows$Percent
  )
  descriptive_rows$Mean_SD_display <- sprintf(
    "%.3f ± %.3f",
    descriptive_rows$Mean,
    descriptive_rows$Standard_deviation
  )
  descriptive_rows$Median_IQR_display <- sprintf(
    "%.3f (%.3f, %.3f)",
    descriptive_rows$Median,
    descriptive_rows$Q1,
    descriptive_rows$Q3
  )

  if (number_groups == 2) {
    # For a two-sample formula call, R labels the statistic W, but its
    # numerical value is the Mann-Whitney U statistic. Do not subtract the
    # rank-sum offset a second time.
    test_result <- wilcox.test(
      Utility ~ Group,
      data = analysis_data,
      exact = FALSE,
      correct = TRUE
    )

    test_name <- "Mann-Whitney U"
    statistic_name <- "U"
    statistic_value <- as.numeric(test_result$statistic)
    test_df <- NA_real_
  } else {
    # Global comparison only; no post-hoc pairwise tests are performed.
    test_result <- kruskal.test(
      Utility ~ Group,
      data = analysis_data
    )

    test_name <- "Kruskal-Wallis"
    statistic_name <- "H"
    statistic_value <- as.numeric(test_result$statistic)
    test_df <- as.numeric(test_result$parameter)
  }

  test_row <- data.frame(
    Variable = display_name,
    Variable_name = variable_name,
    Number_of_groups = number_groups,
    Test = test_name,
    Statistic_name = statistic_name,
    Statistic = statistic_value,
    Degrees_of_freedom = test_df,
    P_value = test_result$p.value,
    Valid_N = valid_n,
    Missing_N = nrow(data) - valid_n,
    stringsAsFactors = FALSE
  )

  test_row$Statistic_display <- sprintf(
    "%.3f",
    test_row$Statistic
  )
  test_row$P_display <- ifelse(
    test_row$P_value < 0.001,
    "<0.001",
    sprintf("%.3f", test_row$P_value)
  )

  list(
    descriptive = descriptive_rows,
    test = test_row
  )
}

univariable_results_list <- lapply(
  analysis_variables,
  function(variable_name) {
    analyse_one_categorical_variable(
      data = dat,
      variable_name = variable_name,
      display_name = unname(variable_labels[variable_name])
    )
  }
)

univariable_descriptive_table <- do.call(
  rbind,
  lapply(univariable_results_list, `[[`, "descriptive")
)

univariable_global_test_table <- do.call(
  rbind,
  lapply(univariable_results_list, `[[`, "test")
)

# Combined table for preparation of manuscript Table 1. The global test is
# repeated across category rows in the CSV so that filtering/sorting does not
# detach a result from its variable.
univariable_descriptive_table$Original_row_order <- seq_len(
  nrow(univariable_descriptive_table)
)

univariable_manuscript_table <- merge(
  univariable_descriptive_table,
  univariable_global_test_table[, c(
    "Variable",
    "Variable_name",
    "Test",
    "Statistic_name",
    "Statistic",
    "Degrees_of_freedom",
    "P_value",
    "Statistic_display",
    "P_display"
  )],
  by = c("Variable_name", "Variable"),
  all.x = TRUE,
  sort = FALSE
)

# Restore the exact variable and category order after merge().
univariable_manuscript_table <-
  univariable_manuscript_table[
    order(univariable_manuscript_table$Original_row_order),
  ]
univariable_manuscript_table$Original_row_order <- NULL
univariable_descriptive_table$Original_row_order <- NULL

# ----------------------------
# 3. Continuous-age exploratory association
# ----------------------------

# Table 1 can retain the prespecified age groups. This additional result checks
# the unadjusted monotonic association using age in its original continuous
# form; it does not replace the continuous age term in the multivariable model.
age_complete <- complete.cases(dat[, c("age", "Utility")])
age_spearman_test <- cor.test(
  dat$age[age_complete],
  dat$Utility[age_complete],
  method = "spearman",
  exact = FALSE
)

continuous_age_association <- data.frame(
  Analysis = "Unadjusted continuous-age association",
  Method = "Spearman rank correlation",
  N = sum(age_complete),
  Spearman_rho = as.numeric(age_spearman_test$estimate),
  P_value = age_spearman_test$p.value,
  P_display = ifelse(
    age_spearman_test$p.value < 0.001,
    "<0.001",
    sprintf("%.3f", age_spearman_test$p.value)
  ),
  stringsAsFactors = FALSE
)

# ----------------------------
# 4. Utility distribution data
# ----------------------------

if (
  min(utility_complete) <= -0.3 ||
    max(utility_complete) > 1
) {
  stop(
    "At least one Utility value lies outside the planned figure range ",
    "(-0.3, 1.0]. Revise utility_breaks before plotting."
  )
}

utility_breaks <- seq(
  from = -0.3,
  to = 1.0,
  by = 0.1
)

format_break <- function(x) {
  ifelse(
    abs(x) < 1e-10,
    "0.0",
    sprintf("%.1f", x)
  )
}

utility_interval_labels <- paste0(
  "(",
  format_break(head(utility_breaks, -1)),
  ",",
  format_break(tail(utility_breaks, -1)),
  "]"
)

utility_interval <- cut(
  utility_complete,
  breaks = utility_breaks,
  labels = utility_interval_labels,
  right = TRUE,
  include.lowest = FALSE
)

if (any(is.na(utility_interval))) {
  stop("At least one Utility value was not assigned to a plotting interval.")
}

utility_frequency <- as.data.frame(
  table(utility_interval),
  stringsAsFactors = FALSE
)
names(utility_frequency) <- c(
  "Utility_interval",
  "Frequency"
)

utility_frequency$Midpoint <- (
  head(utility_breaks, -1) +
    tail(utility_breaks, -1)
) / 2

general_population_mean <- 0.946
maximum_frequency <- max(utility_frequency$Frequency)

# ----------------------------
# 5. Utility distribution figure
# ----------------------------

utility_distribution_plot <- ggplot() +
  geom_col(
    data = utility_frequency,
    aes(
      x = Midpoint,
      y = Frequency
    ),
    width = 0.1,
    fill = "#879BB5",
    color = "white",
    linewidth = 0.4
  ) +
  geom_text(
    data = subset(utility_frequency, Frequency > 0),
    aes(
      x = Midpoint,
      y = Frequency,
      label = Frequency
    ),
    vjust = -0.45,
    size = 3.4,
    family = "Times New Roman"
  ) +
  geom_vline(
    aes(
      xintercept = general_population_mean,
      linetype = "Chinese general population mean"
    ),
    color = "#C00000",
    linewidth = 0.9
  ) +
  # The value is positioned vertically inside the final bar to avoid overlap
  # with the frequency label above that bar.
  annotate(
    geom = "text",
    x = general_population_mean,
    y = maximum_frequency * 0.78,
    label = "0.946",
    color = "#C00000",
    size = 3.4,
    family = "Times New Roman",
    fontface = "bold",
    angle = 90,
    vjust = -0.55
  ) +
  scale_x_continuous(
    breaks = utility_frequency$Midpoint,
    labels = utility_interval_labels,
    limits = c(-0.3, 1.0),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 25),
    expand = expansion(mult = c(0, 0.08))
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "Chinese general population mean" = "dashed"
    )
  ) +
  guides(
    linetype = guide_legend(
      override.aes = list(
        color = "#C00000",
        linewidth = 0.9
      )
    )
  ) +
  labs(
    x = "Health utility values",
    y = "Frequency"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(
      family = "Times New Roman",
      face = "bold",
      size = 10
    ),
    axis.text = element_text(
      family = "Times New Roman",
      size = 9,
      color = "black"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.position = "top",
    legend.justification = "right",
    legend.text = element_text(
      family = "Times New Roman",
      size = 9
    )
  )

print(utility_distribution_plot)

# ================================================================
# 6.Precision of the mean utility estimate using bootstrap
# ================================================================

if (!requireNamespace("boot", quietly = TRUE)) {
  stop(
    "Package 'boot' is required. ",
    "Please install it using install.packages('boot')."
  )
}

# 仅使用非缺失效用值
bootstrap_utility <- dat$Utility[
  !is.na(dat$Utility)
]

# 原始样本均值
observed_mean_utility <- mean(
  bootstrap_utility
)

# bootstrap统计量函数
bootstrap_mean_function <- function(
    data,
    indices
) {
  mean(data[indices])
}

# 固定随机种子，保证结果可以重复
set.seed(20260924)

bootstrap_mean_result <- boot::boot(
  data = bootstrap_utility,
  statistic = bootstrap_mean_function,
  R = 10000
)

bootstrap_mean_result

bootstrap_mean_ci <- boot::boot.ci(
  bootstrap_mean_result,
  conf = 0.95,
  type = c("bca", "perc")
)

bootstrap_mean_ci

# BCa 95%CI
bootstrap_bca_lower <-
  bootstrap_mean_ci$bca[4]

bootstrap_bca_upper <-
  bootstrap_mean_ci$bca[5]

# Percentile 95%CI
bootstrap_percentile_lower <-
  bootstrap_mean_ci$percent[4]

bootstrap_percentile_upper <-
  bootstrap_mean_ci$percent[5]

bootstrap_estimates <- as.numeric(
  bootstrap_mean_result$t
)

bootstrap_standard_error <- sd(
  bootstrap_estimates
)

bootstrap_bias <- mean(
  bootstrap_estimates
) - observed_mean_utility

bootstrap_bca_width <-
  bootstrap_bca_upper -
  bootstrap_bca_lower

bootstrap_bca_half_width <-
  bootstrap_bca_width / 2

bootstrap_precision_table <- data.frame(
  N = length(bootstrap_utility),
  
  Observed_mean = observed_mean_utility,
  
  Bootstrap_resamples = 10000,
  
  Bootstrap_standard_error =
    bootstrap_standard_error,
  
  Bootstrap_bias =
    bootstrap_bias,
  
  BCa_95_CI_lower =
    bootstrap_bca_lower,
  
  BCa_95_CI_upper =
    bootstrap_bca_upper,
  
  BCa_CI_width =
    bootstrap_bca_width,
  
  BCa_CI_half_width =
    bootstrap_bca_half_width,
  
  Percentile_95_CI_lower =
    bootstrap_percentile_lower,
  
  Percentile_95_CI_upper =
    bootstrap_percentile_upper,
  
  stringsAsFactors = FALSE
)

bootstrap_precision_table

bootstrap_precision_table$Mean_BCa_CI_display <- sprintf(
  "%.3f (95%% BCa CI %.3f to %.3f)",
  bootstrap_precision_table$Observed_mean,
  bootstrap_precision_table$BCa_95_CI_lower,
  bootstrap_precision_table$BCa_95_CI_upper
)

bootstrap_precision_table$Mean_BCa_CI_display


# ----------------------------
# 7. Export tables, figure, and caption
# ----------------------------

write.csv(
  overall_utility_summary,
  file.path(output_dir, "01_overall_utility_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  univariable_descriptive_table,
  file.path(output_dir, "02_univariable_descriptive_statistics.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  univariable_global_test_table,
  file.path(output_dir, "03_univariable_global_tests.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  univariable_manuscript_table,
  file.path(output_dir, "04_univariable_manuscript_table.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  continuous_age_association,
  file.path(output_dir, "05_continuous_age_Spearman.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  utility_frequency,
  file.path(output_dir, "06_utility_distribution_frequencies.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  bootstrap_precision_table,
  file.path( output_dir,"07_bootstrap_precision_mean_utility.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

ggsave(
  filename = file.path(output_dir, "Figure3_utility_distribution.png"),
  plot = utility_distribution_plot,
  width = 11,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  filename = file.path(output_dir, "Figure3_utility_distribution.tiff"),
  plot = utility_distribution_plot,
  width = 11,
  height = 6.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)
ggsave(
  filename = file.path(output_dir, "Figure3_utility_distribution.pdf"),
  plot = utility_distribution_plot,
  width = 11,
  height = 6.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

figure_caption <- paste0(
  "Figure 3. Distribution of EQ-5D-5L utility values among participants ",
  "with multidrug-resistant tuberculosis. The dashed red line indicates ",
  "the published mean EQ-5D-5L utility value for the Chinese general ",
  "population. This value is shown for descriptive reference only; no ",
  "inferential comparison or demographic standardisation was performed. ",
  "Utility values below 0 represent health states valued as worse than ",
  "dead under the Chinese EQ-5D-5L value set. One participant had a ",
  "negative utility value."
)

writeLines(
  figure_caption,
  con = file.path(output_dir, "Figure3_caption.txt"),
  useBytes = TRUE
)

message(
  "Univariable analysis and figure completed. Outputs saved to: ",
  normalizePath(output_dir, winslash = "/", mustWork = FALSE)
)
