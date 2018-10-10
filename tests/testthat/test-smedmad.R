### Unit tests for smedmad function

library(CNprep)


### Tests smedmad() results

context("smedmad() results")

test_that("smedmad() must return expected results 01", {
    
    RNGkind("default")
    
    position <- c(1, 6)
    values <- c(0.062073840, 0.119913919, 0.154459489, 0.040994620, -0.082843732,
                0.011907384, 0.001913919, 0.032259489, 0.140994620, -0.000843732)
    
    set.seed(22211)
    
    results <- CNprep:::smedmad(position, values)
    
    expected <- c(0.051534230000000, 0.080065244395500)
    
    expect_equal(results, expected)
})


test_that("smedmad() must return expected results 02", {
    
    RNGkind("default")
    
    position <- c(1, 4)
    values <- c(0.110007384, 0.119913919, 0.854459489, 0.040994620, -0.082843732,
                0.062073840, 0.001913919, 0.032259489, 0.140994620, -0.000843732)
    
    set.seed(3121)
    
    results <- CNprep:::smedmad(position, values)
    
    expected <- c(0.114960651500000, 0.058502876348700)
    
    expect_equal(results, expected)
})

test_that("smedmad() must return expected results when NA present", {
    
    RNGkind("default")
    
    position <- c(1, 5)
    values <- c(0.0520738140, NA, NA, NA, NA,
                0.111907384, NA, NA, NA, -0.000843732)
    
    set.seed(3311)
    
    results <- CNprep:::smedmad(position, values)
    
    expected <- c(0.0520738140000000, 0.000000000000000)
    
    expect_equal(results, expected)
})