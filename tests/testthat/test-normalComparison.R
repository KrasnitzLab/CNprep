### Unit tests for normalComparison.R functions

library(CNprep)


### Tests normalComparison () results

context("normalComparison() results")



test_that("normalComparison() must return expected results 01", {
    
    normalmedian <- c(0.0230326977,  0.0069878681,  0.0013329618,  0.0110395179,  
                      0.0007606011,  0.0023178528, -0.0076653454, -0.0044130592, 
                      -0.0044268312, -0.0362319897)
    
    normalength <- c(11906049, 231977105,  86185990,   2411050, 151410255, 
                     199269234,  34587196,  60357381,  96145497,    620234)
    
    tumormad <- c(0.08502516, 0.13829131, 0.09701612, 0.12886650, 0.09179211, 
                  0.12066915, 0.09360264, 0.09505375, 0.12648441, 0.07690704)
    
    tumormedian <- c(0.046492091, -0.216771966,  0.041802329, -0.227213792,  
                     0.035899255,  0.243933069, -0.032015273, -0.004500177, 
                     -0.284767695, -0.010467050)
    
    tumorerror <- c(0.006267755, 0.009346934, 0.003025827, 0.006431545, 
                    0.003926070, 0.003544793, 0.002222607, 0.003353667, 
                    0.008923441, 0.012611398)
    
    tumorlength <- c(14912327,  21970269,  40526202,  24664080,  38759675, 
                     103077343, 102923669,  40349469,   8897435,   2210770)
    
    RNGkind("default")
    
    set.seed(222)
    
    results <- CNprep:::normalComparison(normalmedian = normalmedian, normalength = normalength, 
                                         tumormedian = tumormedian, 
                                         tumorlength = tumorlength, normalmad = NULL,
                                         normalerror = NULL, tumormad = tumormad,
                                         tumorerror = tumorerror)
                                             
    
    expected <- as.matrix(data.frame(samplesize=c(58, 39, 21, 35, 22, 8, 8, 21, 98, 395),
                           negtail=c(1.0000000000, 0.0000000000, 1.0000000000,
                                     0.0000000000, 1.0000000000, 1.0000000000,
                                     0.0007089442, 0.0402430422, 0.0000000000,
                                     0.0007089442)))
    
    expect_equal(results, expected)
})
