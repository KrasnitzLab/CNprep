### Unit tests for consolidate function

library(CNprep)


data(EMexample)

### Tests consolidate() results

context("consolidate() results")


test_that("consolidate() must return expected results 01", {
    
    RNGkind("default")
    
    set.seed(2111)
    
    expected <- list()
    
    expected[["ngroups"]] <- 4
    expected[["sigmasq"]] <- c(0.0001367895, 0.0001367895,
                               0.0008660118, 0.0001367895)
    expected[["mu"]] <- c(-0.233274515828, -0.081386467651,
                               0.004774756216, 0.245866821040)
    names(expected[["mu"]]) <- c(1:3, 5)
    expected[["pro"]] <- c(0.093264248705, 0.290183786364,
                           0.458520876849, 0.158031088083)
    expected[["groups"]] <- matrix(data = c(1.000000000000, 0.000000000000,
                                            0.000000000000, 0.000000000000,
                                            0.000000000000, 0.000000000000,
                                            1.000000000000, 0.000000000000,
                                            0.000000000000, 0.000000000000,
                                            0.000000000000, 0.000000000000,
                                            1.000000000000, 1.000000000000,
                                            0.000000000000, 0.000000000000,
                                            0.000000000000, 0.000000000000,
                                            0.000000000000, 1.000000000000), 
                                   byrow = TRUE, ncol = 5)
    
    
    results <- CNprep:::consolidate(EMexample, minover = 0.02)
    
    expect_equal(results$ngroups, expected$ngroups)
    expect_equal(results$sigmasq, expected$sigmasq)
    expect_equal(results$mu, expected$mu)
    expect_equal(results$pro, expected$pro)
    expect_equal(results$groups, expected$groups)
})

