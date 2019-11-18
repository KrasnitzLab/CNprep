### Unit tests for applyCNPmask function

library(CNprep)



### Tests applyCNPmask() results

context("applyCNPmask() results")

test_that("applyCNPmask() must return expected results 01", {
    
    segTableTemp01 <- data.frame(ID=c(rep("WZ1", 5)), 
                                 start=c(1, 16, 23, 31, 38),
                                 end=c(15, 22, 30, 37, 50),
                                 num.probes=c(15, 7, 8, 720, 518),
                                 seg.median=c(0.047797239, -0.215466818, 0.043107477, -0.225908644,  
                                              0.037204403),
                                 chrom=c(rep(1, 5)),
                                 chrom.pos.start=c(932544, 16004440, 38093655, 78729960, 103416416),
                                 chrom.pos.end=c(15844870, 37974708, 78619856, 103394039, 142176090),
                                 cytoband.start=c("p36.33", "p36.13", "p34.3", "p31.1", "p21.1"),
                                 cytoband.end=c("p36.13", "p34.3", "p31.1", "p21.1", "q21.1"),
                                 abs.pos.start=c(932544, 16004440, 38093655, 78729960, 103416416),
                                 abs.pos.end=c(15844870, 37974708, 78619856, 103394039, 142176090),
                                 eventIndex=c(0,0,1,0,-1))
    
    cnptableTemp01 <- matrix(c(rep(1, 3), c(932544, 38093688, 123416446), 
                               c(11844870, 48619856, 182176090), rep(1,3)), ncol = 4, byrow = FALSE)
    colnames(cnptableTemp01) <- c("chrom", "start", "end", "cnpindex")
    
    
    RNGkind("default")
    
    set.seed(2111)
    
    results <- applyCNPmask(segTable=segTableTemp01, chrom="chrom", 
                            startPos="chrom.pos.start", endPos="chrom.pos.end",
                            startProbe="start", endProbe = "end", eventIndex="eventIndex", 
                            maskTable = cnptableTemp01, maskChrom="chrom", maskStart="start",
                            maskEnd="end", maskIndex="cnpindex", minCover=0.005,
                            indexVals=c(-1, 1))
    
    expected <- matrix(c(1, 16, 23, 29, 38, 15, 28, 30, 37, 50, 0, 0, 1, 0, 0), ncol = 3, byrow = FALSE)
    colnames(expected) <- c("StartProbe", "EndProbe", "toremove")
    
    expect_equal(results, expected)
})
