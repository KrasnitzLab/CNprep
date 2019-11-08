### Unit tests for breakIntoGaps functions

library(CNprep)



segExample <- data.frame(ID=c(rep("WZ1", 5)), 
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
                         toremove=c(0, 0, 1, 0, 0),
                         stringsAsFactors = FALSE)

segExample2 <- data.frame(ID=c(rep("WZ1", 5)), 
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
                         toremove=c(0, 0, 0, 0, 0),
                         stringsAsFactors = FALSE)

segExample3 <- data.frame(ID=c(rep("WZ1", 5)), 
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
                          toremove=c(0, 0, 0, 0, 1),
                          stringsAsFactors = FALSE)

### Tests breakIntoGaps() results

context("breakIntoGaps() results")

test_that("breakIntoGaps() must return expected results 01", {
  
    
    RNGkind("default")
    
    set.seed(1221)
    
    expected <- matrix(c(1, 16, 23, 28, 38, 15, 27, 30, 37, 50), byrow = FALSE,
                            nrow = 5)
    
    colnames(expected) <- c("start", "end")
   
    results <- CNprep:::breakIntoGaps(segtable = segExample, 
                                      gapind = "toremove",
                                      StartProbe = "start", EndProbe = "end")
    
    expect_equal(results, expected)
      
})


test_that("breakIntoGaps() must return expected results 02", {
    
    
    RNGkind("default")
    
    set.seed(121)
    
    expected <- matrix(c(1, 16, 23, 31, 38, 15, 22, 30, 37, 50), byrow = FALSE,
                       nrow = 5)
    
    colnames(expected) <- c("start", "end")
    
    results <- CNprep:::breakIntoGaps(segtable = segExample2, 
                                      gapind = "toremove",
                                      StartProbe = "start", EndProbe = "end")
    
    expect_equal(results, expected)
    
})



test_that("breakIntoGaps() must return expected results 03", {
    
    
    RNGkind("default")
    
    set.seed(12144)
    
    expected <- matrix(c(1, 16, 23, 31, 38, 15, 22, 30, 50, 50), byrow = FALSE,
                       nrow = 5)
    
    colnames(expected) <- c("start", "end")
    
    results <- CNprep:::breakIntoGaps(segtable = segExample3, 
                                      gapind = "toremove",
                                      StartProbe = "start", EndProbe = "end")
    
    expect_equal(results, expected)
    
})