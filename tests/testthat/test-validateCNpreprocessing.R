### Unit tests for validateCNpreprocessing function

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
                         abs.pos.end=c(15844870, 37974708, 78619856, 103394039, 142176090))


normSegs <- matrix(data = c(0.023032697, 0.0069878681,  0.0013329618, 0.0110395179,  
                            0.0007606011, 0.0023178528, -0.0076653454, -0.0044130592, 
                            -0.0044268312, -0.0362319897, 
                            0.0001533229,  0.0075387269,  0.0078903700,  0.0007193940,  
                            0.0011609480,  0.0001808010, -0.0172001485, -0.0034588337, 
                            -0.0064778796, -0.0001113585), ncol = 1)

colnames(normSegs) <- c("segmedian")

normalLength <- matrix(data = c(11906049, 231977105, 86185990, 2411050, 151410255, 
                                199269234, 34587196, 60357381, 96145497, 620234,
                                179688604, 66855, 35355928, 134792443, 142986004, 
                                15254573, 6675768, 138039242, 38416654, 96961637), ncol=1)

colnames(normalLength) <- c("length")



### Tests validateCNpreprocessing() results

context("validateCNpreprocessing() results")


test_that("validateCNpreprocessing() must return error when nTrial is zero", {
    
    message <- "nTrial must be a positive integer"
    
    expect_error(CNprep:::validateCNpreprocessing(segall=segExample, ratall=NULL, 
                                idCol="ID", startCol="start", endCol="end", 
                                chromCol="chrom", bpStartCol="chrom.pos.start", 
                                bpEndCol="chrom.pos.end", nTrial=0, useEnd=FALSE,
                                blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=1,
                                nJobs=1, modelNames="E", normalLength=normalLength,
                                normalMedian=normSegs), message)
})

test_that("validateCNpreprocessing() must return error when nTrial is negative", {
    
    message <- "nTrial must be a positive integer"
    
    expect_error(CNprep:::validateCNpreprocessing(segall=segExample, ratall=NULL, 
                                                  idCol="ID", startCol="start", endCol="end", 
                                                  chromCol="chrom", bpStartCol="chrom.pos.start", 
                                                  bpEndCol="chrom.pos.end", nTrial=-1, useEnd=FALSE,
                                                  blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=1,
                                                  nJobs=1, modelNames="E", normalLength=normalLength,
                                                  normalMedian=normSegs), message)
})

test_that("validateCNpreprocessing() must return error when nTrial is not an integer", {
    
    message <- "nTrial must be a positive integer"
    
    expect_error(CNprep:::validateCNpreprocessing(segall=segExample, ratall=NULL, 
                                                  idCol="ID", startCol="start", endCol="end", 
                                                  chromCol="chrom", bpStartCol="chrom.pos.start", 
                                                  bpEndCol="chrom.pos.end", nTrial="allo", useEnd=FALSE,
                                                  blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=1,
                                                  nJobs=1, modelNames="E", normalLength=normalLength,
                                                  normalMedian=normSegs), message)
})


test_that("validateCNpreprocessing() must return error when nJobs is zero", {
    
    message <- "nJobs must be a positive integer"
    
    expect_error(CNprep:::validateCNpreprocessing(segall=segExample, ratall=NULL, 
                                                  idCol="ID", startCol="start", endCol="end", 
                                                  chromCol="chrom", bpStartCol="chrom.pos.start", 
                                                  bpEndCol="chrom.pos.end", nTrial=4, useEnd=FALSE,
                                                  blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=1,
                                                  nJobs=0, modelNames="E", normalLength=normalLength,
                                                  normalMedian=normSegs), message)
})


test_that("validateCNpreprocessing() must return error when nJobs is not an integer", {
    
    message <- "nJobs must be a positive integer"
    
    expect_error(CNprep:::validateCNpreprocessing(segall=segExample, ratall=NULL, 
                                                  idCol="ID", startCol="start", endCol="end", 
                                                  chromCol="chrom", bpStartCol="chrom.pos.start", 
                                                  bpEndCol="chrom.pos.end", nTrial=4, useEnd=FALSE,
                                                  blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=1,
                                                  nJobs="NewYork", modelNames="E", normalLength=normalLength,
                                                  normalMedian=normSegs), message)
})

test_that("validateCNpreprocessing() must return error when useEnd is not a logical", {
    
    message <- "useEnd must be a logical value"
    
    expect_error(CNprep:::validateCNpreprocessing(segall=segExample, ratall=NULL, 
                                                  idCol="ID", startCol="start", endCol="end", 
                                                  chromCol="chrom", bpStartCol="chrom.pos.start", 
                                                  bpEndCol="chrom.pos.end", nTrial=4, useEnd="Frida",
                                                  blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=1,
                                                  nJobs=1, modelNames="E", normalLength=normalLength,
                                                  normalMedian=normSegs), message)
})

test_that("validateCNpreprocessing() must return zero when all parameters valid", {
    
    result <- CNprep:::validateCNpreprocessing(segall=segExample, ratall=NULL, 
                                                  idCol="ID", startCol="start", endCol="end", 
                                                  chromCol="chrom", bpStartCol="chrom.pos.start", 
                                                  bpEndCol="chrom.pos.end", nTrial=4, useEnd=FALSE,
                                                  blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=1,
                                                  nJobs=1, modelNames="E", normalLength=normalLength,
                                                  normalMedian=normSegs)
    
    expected <- 0L
    
    expect_equal(result, expected)
})


test_that("validateCNpreprocessing() must return error when bsTimes is not an integer", {
    
    message <- "bsTimes must be a positive integer"
    
    expect_error(CNprep:::validateCNpreprocessing(segall=segExample, ratall=NULL, 
                                                  idCol="ID", startCol="start", endCol="end", 
                                                  chromCol="chrom", bpStartCol="chrom.pos.start", 
                                                  bpEndCol="chrom.pos.end", nTrial=4, useEnd=FALSE,
                                                  blsize=5, minJoin=0.25, cWeight=0.4, bsTimes="hi", chromRange=1,
                                                  nJobs=1, modelNames="E", normalLength=normalLength,
                                                  normalMedian=normSegs), message)
})

test_that("validateCNpreprocessing() must return error when bsTimes is zero", {
    
    message <- "bsTimes must be a positive integer"
    
    expect_error(CNprep:::validateCNpreprocessing(segall=segExample, ratall=NULL, 
                                                  idCol="ID", startCol="start", endCol="end", 
                                                  chromCol="chrom", bpStartCol="chrom.pos.start", 
                                                  bpEndCol="chrom.pos.end", nTrial=4, useEnd=FALSE,
                                                  blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=0, chromRange=1,
                                                  nJobs=1, modelNames="E", normalLength=normalLength,
                                                  normalMedian=normSegs), message)
})

test_that("validateCNpreprocessing() must return error when bsTimes is negative", {
    
    message <- "bsTimes must be a positive integer"
    
    expect_error(CNprep:::validateCNpreprocessing(segall=segExample, ratall=NULL, 
                                                  idCol="ID", startCol="start", endCol="end", 
                                                  chromCol="chrom", bpStartCol="chrom.pos.start", 
                                                  bpEndCol="chrom.pos.end", nTrial=4, useEnd=FALSE,
                                                  blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=-2, chromRange=1,
                                                  nJobs=1, modelNames="E", normalLength=normalLength,
                                                  normalMedian=normSegs), message)
})
