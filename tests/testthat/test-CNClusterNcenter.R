### Unit tests for CNclusterNcenter.R functions

library(CNprep)


data(EMexample)

### Tests CNclusterNcenter() results

context("CNclusterNcenter() results")


test_that("CNclusterNcenter() must return expected results 01", {
    
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
    
    RNGkind("default")
    
    set.seed(11)
    
    profpack <- list()
    profpack[["WZ1"]] <- vector(mode="list", length=4)
    names(profpack[["WZ1"]]) <- c("seg", "rat", "stream", "sub")
    profpack[["WZ1"]]$seg <- segTableTemp01[, c("start", "end", "chrom")]
    dimnames(profpack[["WZ1"]]$seg)[[2]] <- c("StartProbe", 
                                           "EndProbe", "chrom")
    
    profpack[["WZ1"]]$rat <- rep(c(-0.016930740, -0.213708700, 0.001583346,
                               -0.312120554, -0.272878337, 0.193340541, 
                               -0.0001111232, 0.011124328, 0.003334211, 0.0000022222), 5)
    
    profpack[["WZ1"]]$stream <- "WZ1"
    profpack[["WZ1"]]$sub <- 1
    
    
    
    results <- CNprep:::CNclusterNcenter(segrat = profpack[["WZ1"]], blsize = 5, 
                                         minjoin = 0.25, ntrial = 10, bestbic = -1e7, 
                                         modelNames = "E", cweight = 0.4, 
                                         bstimes = 50, chromrange = 1:22)
    
    expected <- matrix(c(-0.0169307400, 0.0000022222, 0.0007927841, -0.0169307400, 0.0000022222, 
                         0.041594444, 0.016489634, 0.009542733, 0.291743003, 0.016489634, 
                         -0.0168762895, 0.0000566727, 0.0008472346, -0.0168762895, 0.0000566727,
                         0.098129264, 0.041742492, 0.038057309, 0.075858316, 0.003258179,
                         0.3800000000, 0.9000000000, 0.9200000000, 0.4400000000, 0.9800000000,
                         0.004777528, 0.031372551, 0.038753772, 0.005801284, 0.038531299,
                         0.6200000, 0.9000000, 0.9200000, 0.5600000, 0.9800000,
                         -0.01772345,  0.000000000, 0.00000000, -0.01772345, 0.000000000,
                         2.833635e-05, 1.467607e-03, 1.467607e-03, 2.833635e-05, 1.467607e-03), 
                        nrow = 5, byrow = FALSE)
    
    colnames(expected) <- c("", "", "mediandev", "segerr", "centerz", 
                            "cpb", "maxz", "maxzmean", "maxzsigma")
    rownames(expected) <- c(1, 2, 2, 1, 2)
    
    expect_equal(results, expected)
})
