### Unit tests for breakIntoCNPs.chrom function

library(CNprep)



### Tests breakIntoCNPs.chrom() results

context("breakIntoCNPs.chrom() results")

test_that("breakIntoCNPs.chrom() must return expected results 01", {
    
    segTableTemp01 <- data.frame(ID=c(rep("WZ1", 5)), 
                                 start=c(1, 15, 23, 31, 37),
                                 end=c(14, 22, 30, 36, 50),
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
                                 eventIndex=c(-1,0,1,0,-1))
    
    cnptableTemp01 <- matrix(c(rep(1, 3), c(932544, 38093688, 123416446), 
                               c(11844870, 48619856, 182176090), rep(1,3)), ncol = 4, byrow = FALSE)
    colnames(cnptableTemp01) <- c("chrom", "start", "end", "cnpindex")
    
    
    RNGkind("default")
    
    set.seed(2116)
    
    results <- CNprep:::breakIntoCNPs.chrom(segtable = segTableTemp01, chrom = "chrom", 
                            startPos = "chrom.pos.start", endPos = "chrom.pos.end",
                            startProbe = "start", endProbe = "end", eventIndex = "eventIndex", 
                            cnptable = cnptableTemp01, cnpchrom = "chrom", cnpstart = "start",
                            cnpend = "end", cnpindex = "cnpindex", mincover = 1,
                            indexvals = c(-1, 1))
    
    expected <- matrix(c(1, 15, 23, 31, 37, 14, 22, 30, 36, 50, 0, 0, 0, 0, 0), ncol = 3, byrow = FALSE)
    colnames(expected) <- c("start", "end", "toremove")
    
    expect_equal(results, expected)
})
