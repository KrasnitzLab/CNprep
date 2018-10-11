### Unit tests for smad function

library(CNprep)


### Tests weighted.median() results

context("weighted.median() results")

test_that("weighted.median() must return expected result 01", {
    
    weights <- c(1, 1, 2, 4, 5)
    values <- c(0.082073840, 0.149913919, 0.144459489, 0.040994620, -0.082843732)
    
    results <- CNprep:::weighted.median(v=values, weights = weights)
    
    expected <- -0.0828437320
    
    expect_equal(results, expected)
})


test_that("weighted.median() must return expected result 02", {
    
    weights <- c(1, 1, 1, 1)
    values <- c(0.082073840, 0.149913919, 0.144459489, 0.040994620)
    
    results <- CNprep:::weighted.median(v=values, weights = weights)
    
    expected <- 0.082073840
    
    expect_equal(results, expected)
})