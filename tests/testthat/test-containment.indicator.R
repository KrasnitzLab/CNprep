### Unit tests for containment.indicator.R functions

library(CNprep)


### Tests containment.indicator() results

context("containment.indicator() results")

test_that("containment.indicator() must return expected results 01", {
    
    ustart <- 12842827
    uend <- 13456522
    dstart <- ustart
    dend <- uend
    
    
    expected <- matrix(data=c(1, 1), nrow=1, 
                       dimnames=list(c(NULL), c("startafterstart", "endbeforeend")))
    
    results <- CNprep:::containment.indicator(vstart=ustart, vend=uend, 
                                        wstart=dstart, wend=dend)
    
    expect_equal(results, expected)
})
