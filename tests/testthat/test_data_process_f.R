library(testthat)
library(SurvPredPipe)  # Load your package
# Load the function
source("/Users/kaurh8/Documents/Survival_Pred_Package/Surv_Pred_Functions/SurvPredPipe/R/data_process_f.R")

# Define test cases
test_that("data_process_f function works as expected", {
  
  # Test case 1: Ensure function works with valid input
  data <- "./extdata/TCGA-LGG_protein_coding_FPKM_data_with_clin_data.txt"
  col_num <- 20
  surv_time <- "OS.time"
  output <- "New_data.txt"
  
  # Run the function
  data_process_f(data, col_num, surv_time, output)
  
  # Assert that the output file was created
  expect_true(file.exists(output))
  
  # Cleanup: remove the output file
  file.remove(output)
  
  # Test case 2: Ensure function raises an error with empty input
  expect_error(data_process_f(data = NULL, col_num = NULL, surv_time = NULL, output = NULL),
               "Error: Empty input variable detected.")
  
  # Test case 3: Ensure function raises an error with missing values in input
  expect_error(data_process_f(data = data, col_num = col_num, surv_time = NULL, output = output),
               "Error: Missing values in input variables.")
})

# Run the tests
test_file("test_data_process_f.R")
