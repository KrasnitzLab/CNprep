### Unit tests for breakIntoGaps function

library(CNprep)



### Tests breakIntoGaps() results

context("breakIntoGaps() results")

test_that("breakIntoGaps() must return expected results 01", {
    
    RNGkind("default")
    
    set.seed(211231)
    
    segTableTemp01 <- data.frame(ID=c(rep("WZ1", 5)), 
                                 start=c(1, 13, 23, 31, 38),
                                 end=c(12, 21, 30, 37, 50),
                                 num.probes=c(15, 7, 8, 720, 518),
                                 seg.median=c(0.047788239, -0.215444818, 0.043107477, -0.225908644,  
                                              0.037204403),
                                 chrom=c(rep(1, 5)),
                                 chrom.pos.start=c(932544, 16004440, 38093655, 78729960, 103416416),
                                 chrom.pos.end=c(15844870, 37974706, 78619856, 103394039, 142176090),
                                 cytoband.start=c("p36.33", "p36.13", "p34.3", "p31.1", "p21.1"),
                                 cytoband.end=c("p36.13", "p34.3", "p31.1", "p21.1", "q21.1"),
                                 abs.pos.start=c(932544, 16004440, 38093655, 78729960, 103416416),
                                 abs.pos.end=c(15844870, 37974706, 78619856, 103394039, 142176090),
                                 eventIndex=c(0,0,1,0,-1),
                                 toremove=c(0,0,1,0,0))
    
    action <- "toremove"
    startID <- "start"
    endID <- "end"
    
    results <- CNprep:::breakIntoGaps(segtable =  segTableTemp01, gapind = action, 
                                      StartProbe = startID, EndProbe = endID)
    
    expected <- matrix(c(1, 13, 23, 26, 38, 12, 25, 30, 37, 50), ncol=2, byrow = FALSE)
    colnames(expected) <- c("start", "end")
    
    expect_equal(results, expected)
})

