library(testthat)
library(SurvPredPipe)  # Load your package

# Load the function
source("/Users/kaurh8/Documents/Survival_Pred_Package/Surv_Pred_Functions/SurvPredPipe/R/Nomogram_generate_f.R")
# Define the test cases
test_that("Nomogram_generate_f function test", {
  # Test case 1: Check for empty input variables
  expect_error(Nomogram_generate_f(
    data = NULL,
    Feature_List = "feature_list_for_Nomogram.txt",
    surv_time = "OS_month",
    surv_event = "OS"
  ), "Error: Empty input variable detected.")

  # Test case 2: Check for missing values in input variables
  expect_error(Nomogram_generate_f(
    data = "Train_Data_Nomogram_input.txt",
    Feature_List = NA,
    surv_time = "OS_month",
    surv_event = "OS"
  ), "Error: Missing values in input variables.")

  # Add more test cases as needed
})


test_file("test_Nomogram_generate_f.R")
