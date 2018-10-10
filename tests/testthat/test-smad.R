### Unit tests for smad function

library(CNprep)


### Tests smad() results

context("smad() results")

test_that("smad() must return expected results 01", {
    
    RNGkind("default")
    
    position <- c(1, 5)
    values <- c(0.082073840, 0.149913919, 0.144459489, 0.040994620, -0.082843732,
                0.111907384, 0.001913919, 0.032259489, 0.140994620, -0.000843732)
    
    set.seed(2111)
    
    results <- CNprep:::smad(position, values)
    
    expected <- 0.092492963207400
    
    expect_equal(results, expected)
})


test_that("smad() must return expected results 02", {
    
    RNGkind("default")
    
    position <- c(1, 5)
    values <- c(0.072073840, 0.219913919, 0.154459489, 0.040994620, -0.082843732,
                0.111907384, 0.001913919, 0.032259489, 0.140994620, -0.000843732)
    
    set.seed(32211)
    
    results <- CNprep:::smad(position, values)
    
    expected <- 0.122144963207400
    
    expect_equal(results, expected)
})

test_that("smad() must return expected results when NA present", {
    
    RNGkind("default")
    
    position <- c(1, 5)
    values <- c(0.07211073840, NA, NA, NA, NA,
                0.111907384, NA, NA, NA, -0.000843732)
    
    set.seed(33311)
    
    results <- CNprep:::smad(position, values)
    
    expected <- 0.000000000000000
    
    expect_equal(results, expected)
})