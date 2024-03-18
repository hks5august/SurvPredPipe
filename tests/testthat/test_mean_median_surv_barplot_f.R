library(testthat)
library(SurvPredPipe)  # Load your package


# Load the function
source("/Users/kaurh8/Documents/Survival_Pred_Package/Surv_Pred_Functions/SurvPredPipe/R/mean_median_surv_barplot_f.R")

# Define the test cases
test_that("mean_median_surv_barplot_f function test", {
  # Test case 1: Check for empty input variables
  expect_error(mean_median_surv_barplot_f(
    surv_mean_med_data = NULL,
    selected_sample = "TCGA-DH-5140-01"
  ), "Error: Empty input variable detected.")

  # Test case 2: Check for missing values in input variables
  expect_error(mean_median_surv_barplot_f(
    surv_mean_med_data = "mean_median_survival_time_data.txt",
    selected_sample = NA
  ), "Error: Missing values in input variables.")

  # Add more test cases as needed
})

test_file("test_mean_median_surv_barplot_f.R")

