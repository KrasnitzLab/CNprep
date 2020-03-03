### Unit tests for get.center function

library(CNprep)



### Tests get.center() results

context("get.center() results")

test_that("get.center() must return expected results 01", {
    
    demoEM <- list()
    demoEM[["mu"]] <- c(-0.23626, -0.08108, -0.02205, 0.03059, 0.24482)
    demoEM[["pro"]] <- rep(0.2, 5)
    demoEM[["z"]] <- matrix(data=c(1.19e-118, 2.81e-25, 5.87e-08, 9.99e-1,  
         1.86e-52, 2.03e-117, 9.19e-25, 1.02e-07, 9.99e-01, 1.92e-53, 1.00e+0, 
         1.34e-23, 1.72e-50, 1.08e-82, 6.45e-295, 1.00e+00, 1.39e-20, 2.51e-46, 
         1.67e-77, 1.47e-285, 8.86e-63, 1.21e-04, 9.99e-01, 1.89e-05, 7.93e-106,
         7.59e-60, 7.76e-04, 9.99e-01, 3.60e-06, 1.75e-109, 0.00e+0, 1.61e-147, 
         1.08e-98, 2.31e-63, 1.00e+0, 0.00e+0, 1.18e-147, 8.37e-99, 1.88e-63, 
         1.00e+0, 3.51e-75, 9.79e-01, 4.55e-08, 2.06-02, 2.14e-90, 7.07e-79,
         8.58e-01, 3.96e-09,  1.41e-01, 6.42e-86), ncol=5, byrow=TRUE)
    demoEM[["groups"]] <- diag(x=1, nrow=5, ncol=5, names=TRUE)
    demoEM[["ngroups"]] <- 5
    demoEM[["sigmasq"]] <- rep(1.533e-3, 5)
    
    results <- CNprep:::get.center(emfit=demoEM, minCenter=0.30)
    
    expected <- list()
    expected[["mu"]] <- c(-0.23626, -0.08108, 0.00427, 0.24482)
    expected[["pro"]] <- c(0.2, 0.2, 0.4, 0.2)
    expected[["z"]] <- matrix(data=c(1.19e-118, 2.81e-25, 9.990000587e-1, 1.86e-52, 
                                     2.03e-117, 9.19e-25, 9.99000102e-1, 1.92e-53, 
                                     1.00e+0,   1.34e-23, 1.72e-50, 6.45e-295, 
                                     1.00e+00, 1.39e-20, 2.51e-46, 1.47e-285, 
                                     8.86e-63, 1.21e-04, 9.990189e-1, 7.93e-106,
                                     7.59e-60, 7.76e-04, 9.990036e-1, 1.75e-109, 
                                     0.00e+0, 1.61e-147, 2.31e-63, 1.00e+0, 
                                     0.00e+0, 1.18e-147, 1.88e-63, 1.00e+0, 
                                     3.51e-75, 9.79e-1, 6.000004550e-2, 2.14e-90, 
                                     7.07e-79, 8.58e-1, 1.4100000396e-1, 6.42e-86), ncol=4, byrow=TRUE)
    expected[["groups"]] <- diag(x=1, nrow=4, ncol=5, names=TRUE)
    expected[["groups"]][3,4] <- 1
    expected[["groups"]][4,4] <- 0
    expected[["groups"]][4,5] <- 1
    expected[["ngroups"]] <- 4
    expected[["sigmasq"]] <- c(1.533e-3, 1.533e-3, 2.2257424e-3, 1.533e-3)
    expected[["center"]] <- 3
    
    expect_equal(results, expected)
})
