### Unit tests for segsample.R functions

library(CNprep)


### Tests segsample() results

context("segsample() results")

test_that("segsample() must return expected results", {
    
    segData <- data.frame(StartProbe=c(1, 7, 13),
                          EndProbe=c(6, 12, 15),
                          chrom=c(1,1,1),
                          segmedian=c(0.08662475, 0.07319237, 0.07689544),
                          segmad=c(0.06404208, 0.04764233, 0.08833202))
    ratcol <- c(0.072073840, 0.119913919,  0.154459489,  0.040994620, -0.082843732, 
                0.093052725,  0.170908930,  0.100289490,  0.086624752, -0.003855011, 
                -0.195791649, 0.063634112,  0.109449474,  0.043428961,  0.160174529)
    
    
    RNGkind("default")
    
    set.seed(222)
    
    results <- CNprep:::segsample(mysegs=segData, ratcol=ratcol, startcol="StartProbe", endcol="EndProbe",
                                 blocksize=5,times=0)
    
    colnames(results)[3] <- "results"
    
    expected <- data.frame(StartProbe=c(1, 7),
                           EndProbe=c(6, 12),
                           results=c(0.09305272500, 0.08662475200))
    
    expect_equal(results, expected)
})