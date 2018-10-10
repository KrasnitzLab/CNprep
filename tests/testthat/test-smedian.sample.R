### Unit tests for smedian.sample function

library(CNprep)


### Tests CNpreprocessing() results

context("smedian.sample() results")

test_that("smedian.sample() must return expected results 01", {
    
    RNGkind("default")
    
    position <- c(1, 5)
    values <- c(0.072073840, 0.119913919, 0.154459489, 0.040994620, -0.082843732,
                0.111907384, 0.001913919, 0.032259489, 0.140994620, -0.000843732)
    
    set.seed(2111)
    
    results <- CNprep:::smedian.sample(position, values)
    
    expected <- 0.0720738400
    
    expect_equal(results, expected)
})


test_that("smedian.sample() must return expected results 02", {
    
    RNGkind("default")
    
    position <- c(1, 5)
    values <- c(0.072073840, 0.119913919, 0.154459489, 0.040994620, -0.082843732,
                0.111907384, 0.001913919, 0.032259489, 0.140994620, -0.000843732)
    
    set.seed(311)
    
    results <- CNprep:::smedian.sample(position, values)
    
    expected <- 0.0409946200
    
    expect_equal(results, expected)
})