### Unit tests for validateMakeCNPmask function

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



### Tests validateMakeCNPmask() results

context("validateMakeCNPmask() results")


test_that("validateMakeCNPmask() must return error when nProf is zero", {
    
    message <- "nProf must be a positive integer"
    
    expect_error(CNprep:::validateMakeCNPmask(imat=segExample, chromCol="chrom",
                                startCol="start", endCol="end", nProf=0, 
                                uThresh=0.20, dThresh=0.10), message)
})

test_that("validateMakeCNPmask() must return error when dThresh is a string", {
    
    message <- "dThresh must be a ratio between 0 and 1"
    
    expect_error(CNprep:::validateMakeCNPmask(imat=segExample, chromCol="chrom",
                                              startCol="start", endCol="end", nProf=10, 
                                              uThresh=0.10, dThresh="Chicago"), message)
})

test_that("validateMakeCNPmask() must return error when uThresh is a string", {
    
    message <- "uThresh must be a ratio between 0 and 1"
    
    expect_error(CNprep:::validateMakeCNPmask(imat=segExample, chromCol="chrom",
                                              startCol="start", endCol="end", nProf=10, 
                                              uThresh="Cleveland", dThresh=0.4), message)
})

test_that("validateMakeCNPmask() must return error when dThresh negative", {
    
    message <- "dThresh must be between 0 and 1"
    
    expect_error(CNprep:::validateMakeCNPmask(imat=segExample, chromCol="chrom",
                                              startCol="start", endCol="end", nProf=10, 
                                              uThresh=0.10, dThresh=-0.01), message)
})


test_that("validateMakeCNPmask() must return error when uThresh negative", {
    
    message <- "uThresh must be between 0 and 1"
    
    expect_error(CNprep:::validateMakeCNPmask(imat=segExample, chromCol="chrom",
                                              startCol="start", endCol="end", nProf=10, 
                                              uThresh=-0.01, dThresh=0.01), message)
})


test_that("validateMakeCNPmask() must return error when dThresh superior to uThresh", {
    
    message <- "uThresh must be equal or superior to dThresh"
    
    expect_error(CNprep:::validateMakeCNPmask(imat=segExample, chromCol="chrom",
                                              startCol="start", endCol="end", nProf=10, 
                                              uThresh=0.10, dThresh=0.20), message)
})


test_that("validateMakeCNPmask() must return zero when all parameters valid", {
    
    expected <- 0L
    
    result <- CNprep:::validateMakeCNPmask(imat=segExample, chromCol="chrom",
                                startCol="start", endCol="end", nProf=10, 
                                uThresh=0.20, dThresh=0.10)
    
    expect_equal(result, expected)
})
