### Unit tests for smedian function

library(CNprep)


### Tests smedian() results

context("smedian() results")

test_that("smedian() must return expected results 01", {
    
    RNGkind("default")
    
    position <- c(1, 5)
    values <- c(0.062073840, 0.119913919, 0.154459489, 0.040994620, -0.082843732,
                0.111907384, 0.001913919, 0.032259489, 0.140994620, -0.000843732)
    
    set.seed(2211)
    
    results <- CNprep:::smedian(position, values)
    
    expected <- 0.0620738400
    
    expect_equal(results, expected)
})


test_that("smedian() must return expected results 02", {
    
    RNGkind("default")
    
    position <- c(1, 4)
    values <- c(0.110007384, 0.119913919, 0.854459489, 0.040994620, -0.082843732,
                0.062073840, 0.001913919, 0.032259489, 0.140994620, -0.000843732)
    
    set.seed(311)
    
    results <- CNprep:::smedian(position, values)
    
    expected <- 0.114960651500
    
    expect_equal(results, expected)
})

test_that("smedian() must return expected results when NA present", {
    
    RNGkind("default")
    
    position <- c(1, 5)
    values <- c(0.052073840, NA, NA, NA, NA,
                0.111907384, NA, NA, NA, -0.000843732)
    
    set.seed(3311)
    
    results <- CNprep:::smedian(position, values)
    
    expected <- 0.05207384000
    
    expect_equal(results, expected)
})