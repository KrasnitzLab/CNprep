### Unit tests for breakIntoCNPs.chrom.R functions

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

rateExample <- matrix(data = c(0.0720738404241163, 0.119913919265996, 0.154459489283567, 0.0409946196196852, -0.0828437318817242, 
                               0.0930527245776551, 0.170908929806745, 0.100289489712059, 0.08662475225992, -0.00385501101587094, -0.195791648504205, 
                               0.0636341122831496, 0.109449474459845, 0.043428961427158, 0.160174528536138, 0.0410580626539806, 0.0405224616647681, 
                               0.0986143513758616, 0.0731923745507058, 0.0121261797978406, 0.107795728119586, 0.0832422912347904, 0.0727354592601023, 
                               0.1475937155592, 0.0810554254709479, 0.162564668451733, 0.0590673628590426, -0.0728672744452223, 0.0284354530744622, 
                               0.160126815673499, 0.190799021333644, -0.0269711432208357, 0.0203580612494725, -0.00731229981882388, -0.0440263312613861, 
                               -0.0817496261303279, 0.0146393203392997, 0.0860431917828892, 0.0213604328945665, 0.0789323375946298, 0.0757082196458755, 
                               -0.0507940506094144, 0.0217577735749916, 0.0281455710909498, 0.0196702655955798, 0.154422985605817, 0.260452283950961, 0.122065762785514, 
                               -0.0197506768456278, -0.0277870900508885), ncol = 1)

colnames(rateExample) <- c("WZ1")

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




### Tests breakIntoCNPs.chrom() results

context("breakIntoCNPs.chrom() results")

test_that("breakIntoCNPs.chrom() must return expected results 01", {
    
    
    segTableTemp01 <- data.frame(ID=c(rep("WZ1", 5)), 
                                 start=c(1, 16, 23, 31, 38),
                                 end=c(15, 22, 30, 37, 50),
                                 num.probes=c(15, 7, 8, 720, 518),
                                 seg.median=c(0.037797239, -0.515466818, 0.143107477, -0.225908644,  
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
                               c(11844870, 48619856, 182176090), rep(1,3)), ncol=4, byrow=FALSE)
    
    colnames(cnptableTemp01) <- c("chrom", "start", "end", "cnpindex")
    
    
    RNGkind("default")

    set.seed(111223)
    
    results <- CNprep:::breakIntoCNPs.chrom(segTable=segTableTemp01, chrom="chrom", 
                            startPos="chrom.pos.start", endPos="chrom.pos.end", 
                            startProbe="start",
                            endProbe="end", 
                            eventIndex="eventIndex", 
                            cnpTable=cnptableTemp01, cnpChrom="chrom", 
                            cnpStart="start", cnpEnd="end", cnpIndex="cnpindex", 
                            minCover=0.005, indexVals=c(-1, 1))
    
    
    expected <- matrix(c(1, 16, 23, 25, 38, 15, 24, 30, 37, 50, 0, 0, 1, 0, 0), ncol=3, byrow=FALSE)
    colnames(expected) <- c("start", "end", "toremove")
    
    expect_equal(results, expected)
    
})


test_that("breakIntoCNPs.chrom() must return expected results when eventIndex all zero", {
    
    
    segTableTemp01 <- data.frame(ID=c(rep("WZ1", 5)), 
                                 start=c(1, 16, 23, 31, 38),
                                 end=c(15, 22, 30, 37, 50),
                                 num.probes=c(15, 7, 8, 720, 518),
                                 seg.median=c(0.037797239, -0.515466818, 0.143107477, -0.225908644,  
                                              0.037204403),
                                 chrom=c(rep(2, 5)),
                                 chrom.pos.start=c(932544, 16004440, 38093655, 78729960, 103416416),
                                 chrom.pos.end=c(15844870, 37974708, 78619856, 103394039, 142176090),
                                 cytoband.start=c("p36.33", "p36.13", "p34.3", "p31.1", "p21.1"),
                                 cytoband.end=c("p36.13", "p34.3", "p31.1", "p21.1", "q21.1"),
                                 abs.pos.start=c(932544, 16004440, 38093655, 78729960, 103416416),
                                 abs.pos.end=c(15844870, 37974708, 78619856, 103394039, 142176090),
                                 eventIndex=c(0,0,0,0,0))
    
    cnptableTemp01 <- matrix(c(rep(1, 3), c(932544, 38093688, 123416446), 
                               c(11844870, 48619856, 182176090), rep(1,3)), ncol=4, byrow=FALSE)
    
    colnames(cnptableTemp01) <- c("chrom", "start", "end", "cnpindex")
    
    
    RNGkind("default")
    
    set.seed(123)
    
    results <- CNprep:::breakIntoCNPs.chrom(segTable=segTableTemp01, chrom="chrom", 
                                            startPos="chrom.pos.start", endPos="chrom.pos.end", 
                                            startProbe="start",
                                            endProbe="end", 
                                            eventIndex="eventIndex", 
                                            cnpTable=cnptableTemp01, cnpChrom="chrom", 
                                            cnpStart="start", cnpEnd="end", cnpIndex="cnpindex", 
                                            minCover=0.005, indexVals=c(-1, 1))
    
    
    expected <- matrix(c(1, 16, 23, 31, 38, 15, 22, 30, 37, 50, 0, 0, 0, 0, 0), ncol=3, byrow=FALSE)
    colnames(expected) <- c("start", "end", "toremove")
    
    expect_equal(results, expected)
    
})
