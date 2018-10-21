### Unit tests for normalComparison.R functions

library(CNprep)


### Tests normalComparison() results

context("normalComparison() results")

test_that("normalComparison() must return expected results 01", {
    
    RNGkind("default")
    
    set.seed(112211)
    
    normalmedian <- matrix(c(0.0230326970, 0.0069878681, 0.0013329618, 
                                0.0110395179, 0.0007606011), ncol=1)
    normalength <- matrix(c(11906049, 231977105, 86185990, 2411050, 
                                151410255), ncol=1)
    
    tumormedian <- c(0.066954487, 0.053522109, 0.057225177, -0.026982565, 
                        0.008475305)
    tumorlength <- c(14912327, 21970269, 40526202, 24664080, 38759675)
    
    tumormad <- c(0.06404208, 0.04764233, 0.08833202, 0.04102408, 0.07529646)
    tumorerror <- c(NA, NA, NA, NA, NA)
    
    results <- CNprep:::normalComparison(normalmedian=normalmedian, 
                                normalength=normalength, tumormedian=tumormedian,
                                tumorlength=tumorlength, normalmad=NULL,
                                normalerror=NULL, tumormad=tumormad, 
                                tumorerror=tumorerror)

    
    expected <- matrix(c(32, 22, 11, 19, 12, 1.000000000000, 1.0000000000000, 
                         1.000000000000, 0.0000000000000, 0.970412519962757), 
                         byrow = FALSE, ncol=2)
    colnames(expected) <- c("samplesize", "negtail")
    
    expect_equal(results, expected)
})