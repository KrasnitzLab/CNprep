### Unit tests for makeCNPmask function

library(CNprep)



### Tests makeCNPmask() results

context("makeCNPmask() results")

test_that("makeCNPmask() must return expected results 01", {
    
    ampTemp01 <- data.frame(copy.num=c(rep("amp", 6)), 
                                 chrom=c(1, 2,  1,  5,  1, 3),
                                 chrom.start=c(12842827, 86198238, 13071912, 233492, 13071912, 199030879),
                                 chrom.end=c(13439275, 86575258, 13439275, 292937, 13456522, 199326099), stringsAsFactors = FALSE)
    
    RNGkind("default")
    
    set.seed(211)
    
    results <- makeCNPmask(imat=ampTemp01, chromCol="chrom", 
                           startCol="chrom.start", endCol="chrom.end",
                           nProf=3, uthresh=0.02, dthresh=0.008)
    
    expected <- matrix(c(1, 2, 3, 5, 
                         12842827, 86198238, 199030879, 233492, 
                         13456522, 86575258, 199326099, 292937), 
                       ncol = 3, byrow = FALSE)
    colnames(expected) <- c("chrom", "start", "end")
    
    expect_equal(results, expected)
})


test_that("makeCNPmask() must return no result", {
    
    ampTemp01 <- data.frame(copy.num=c(rep("amp", 6)), 
                            chrom=c(1, 2,  4,  5,  6, 8),
                            chrom.start=c(12842827, 86198238, 13071912, 233492, 13071912, 199030879),
                            chrom.end=c(13439275, 86575258, 13439275, 292937, 13456522, 199326099), stringsAsFactors = FALSE)
    
    RNGkind("default")
    
    set.seed(2211)
    
    results <- makeCNPmask(imat=ampTemp01, chromCol="chrom", 
                           startCol="chrom.start", endCol="chrom.end",
                           nProf=100, uthresh=0.02, dthresh=0.008)
    
    expected <- matrix(data = as.double(c(NULL)), ncol = 3, nrow = 0,
                       dimnames = list(NULL, c("chrom", "start", "end")))
    
    expect_equal(results, expected)
})

