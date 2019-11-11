### Unit tests for getz.R functions

library(CNprep)


data(EMexample)

### Tests getz() results

context("getz() results")

test_that("getz() must return expected results 01", {
    
    
    RNGkind("default")
    
    set.seed(2411)
    
    logr <- c(0.09305272,  0.04052246, 0.16012682, -0.04402633, 0.07570822)
    
    zgroup <- matrix(data = rep(0, 4 * 5), nrow = 4)
    zgroup[1, 1:2] <- 1
    zgroup[2, 3] <- 1
    zgroup[3, 4] <- 1
    zgroup[4, 5] <- 1
    
    
    
    expected <- matrix(data = c(4.7368048124e-43, 4.2288297443e-24, 0.0000000000e+00, 
                                4.4224316883e-02, 0.0000000000e+00, 8.9968641384e-16,
                                9.1795728352e-07, 0.0000000000e+00, 9.5577567982e-01,
                                8.4843243542e-13, 1.0000000000e+00, 9.9999908204e-01,
                                7.1054273576e-15, 3.2942020223e-09, 1.0000000000e+00,
                                4.4237131507e-32, 0.0000000000e+00, 1.0000000000e+00,
                                0.0000000000e+00, 0.0000000000e+00), byrow = FALSE, nrow=5)
    
    
    results <- CNprep:::getz(logr = logr, emfit = EMexample, zgroup = zgroup, 
                                times = 1)
    
    
    expect_equal(results, expected)
})


test_that("getz() must return expected results 02", {
    
    
    RNGkind("default")
    
    set.seed(444)
    
    tempEM <- EMexample
    tempEM$z <- NULL
    tempEM$parameters$mean <- tempEM$parameters$mean[1]
    tempEM$parameters$pro <- tempEM$parameters$pro[1]
    
    
    logr <- c(0.01305272,  0.01052246, 0.18012682, -0.09402633, 0.07570822)
    
    
    zgroup <- matrix(data = rep(0, 1 * 5), nrow = 5)
    zgroup[1, 1] <- 1 
    
    
    expected <- matrix(data = c(rep(1.0000000000e+00, 5), 
                                rep(0.0000000000e+00, 20)), byrow = FALSE, nrow=5)
    
    
    results <- CNprep:::getz(logr = logr, emfit = tempEM, zgroup = zgroup, 
                             times = 1)
    
    
    expect_equal(results, expected)
})


