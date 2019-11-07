### Unit tests for centerprob.R functions

library(CNprep)

data("EMexample")

### Tests centerprob() results

context("centerprob() results")

test_that("centerprob() must return expected results 01", {
    
    logr <- c(0.09305272,  0.04052246, 0.16012682, -0.04402633, 0.07570822)
    
    zgroup <- matrix(data = rep(0, 4 * 5), nrow = 4)
    zgroup[1, 1:2] <- 1
    zgroup[2, 3] <- 1
    zgroup[3, 4] <- 1
    zgroup[4, 5] <- 1
    
    results <- CNprep:::centerprob(logr = logr, emfit = EMexample, 
                                   zgroup = zgroup, times = 1, center =  1)
    
    expected <- c(9.9908332267e-51, 7.337291e-26, 0.000000e+00, 
                  5.303009e-04, 0.000000e+00)
    
    expect_equal(results, expected)
})

test_that("centerprob() must return expected results 02", {
    
    
    EmTemp <- EMexample
    EmTemp$parameters$mean <- EmTemp$parameters$mean[1]
    EmTemp$parameters$pro <- EmTemp$parameters$pro[1]
    
    
    logr <- c(0.01305272,  0.00052246, 0.06012682, -0.04102633, 0.01510822)
    
    zgroup <- matrix(data = rep(0, 1 * 5), nrow = 5)
    zgroup[1, 1] <- 1 
    
    results <- CNprep:::centerprob(logr = logr, emfit = EmTemp, 
                                   zgroup = zgroup, times = 1, center =  1)
    
    expected <- c(9.0002107915e-99, 3.3657353567e-89, 0.000000e+00, 
                  5.1539497408e-61, 0.000000e+00)
    
    expect_equal(results, expected)
})