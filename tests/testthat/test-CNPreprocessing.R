### Unit tests for CNVMetricsMethods.R functions

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



### Tests CNpreprocessing() results

context("CNpreprocessing() results")

test_that("CNpreprocessing() must return expected results 01", {
    
    set.seed(444)
    
    
    results <- CNpreprocessing(segall=segExample, ratall=rateExample, idCol="ID", 
        "start", "end", chromCol="chrom", bpStartCol="chrom.pos.start", 
        bpEndCol="chrom.pos.end",
        blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=1, nJobs=1,
        modelNames="E", normalLength=normalLength,
        normalMedian=normSegs)
    row.names(results) <- NULL
    
    expected <- segExample
    
    expected$segmedian <- c(0.0866247523, 0.0731923746, 0.0768954424, -0.0073122998, 0.0281455711)
    expected$segmad <- c(0.0640420795, 0.0476423308, 0.0883320200, 0.0410240773, 0.0752964600)
    expected$mediandev <- c(0.022990639976770, 0.009558262267556, 
                            0.013261330082376, -0.070946412101973, 
                            -0.035488541192200)
    expected$segerr <- as.double(rep(NA, 5))
    expected$centerz <- c(1, 1, 1, 1, 1)
    expected$marginalprob <- c(0.000000000000000, 0.000000000000000, 
                               0.000139724182812, 0.000000000000000, 
                               0.000000000000000)
    expected$maxz <- rep(1, 5)
    expected$maxzmean <- c(0.000000000000000, 0.000000000000000, 
                           0.000000000000000, 0.000000000000000, 
                           0.000000000000000)
    expected$maxzsigma <- c(0.045757298969319, 0.045757298969319, 
                            0.045757298969319, 0.045757298969319, 
                            0.045757298969319)
    expected$samplesize <- c(111, 75, 41, 67, 42)
    expected$negtail <- c(0.992841083588778, 0.991391357872609, 
                          0.992841083588778, 0.000000000000000, 
                          0.000372936761926)
    row.names(expected) <- NULL
    
    expect_equal(results, expected)
})

test_that("CNpreprocessing() must return expected results 02", {
    
    set.seed(424)
    
    results <- CNpreprocessing(segall=segExample, ratall=rateExample, idCol="ID", 
                                "start", "end", chromCol="chrom", 
                                bpStartCol="chrom.pos.start", 
                                bpEndCol="chrom.pos.end",
                                blsize=3, minJoin=0.25, cWeight=0.4, bsTimes=1, 
                                chromRange=1, modelNames="E", 
                                normalLength=normalLength,
                                normalMedian=normSegs)
    row.names(results) <- NULL
    
    expected <- segExample
    
    expected$segmedian <- c(0.0866247523, 0.0731923746, 0.0768954424, 
                            -0.0073122998, 0.0281455711)
    expected$segmad <- c(0.064042079488653, 0.047642330818085, 
                         0.088332019979936, 0.041024077319856, 
                         0.075296460018356)
    expected$mediandev <- c(0.064866978684928, 0.051434600975714, 
                            0.055137668790534, -0.029070073393815, 
                            0.006387797515958)
    expected$segerr <- as.double(rep(NA, 5))
    expected$centerz <- c(1, 0, 0, 1, 0)
    expected$marginalprob <- c(0.109372964234594, 0.000000000000000, 
                               0.000000000000000, 0.008601807590315, 
                               0.000000000000000)
    expected$maxz <- c(1.000000000000000, 1.000000000000000, 
                       0.538541655965071, 1.000000000000000, 
                       0.990162113892837)
    expected$maxzmean <- c(0.000000000000000, 0.059475806154862, 
                           0.048267482678576, 0.000000000000000, 
                           0.048267482678576)
    expected$maxzsigma <- c(0.016795531137887, 0.002310603657899, 
                            0.002310603657899, 0.016795531137887, 
                            0.002310603657899)
    expected$samplesize <- c(111, 75, 41, 67, 42)
    expected$negtail <- c(1.000000000000000, 1.000000000000000, 
                          1.000000000000000, 0.000372936761926, 
                          0.830608091605016)
    row.names(expected) <- NULL
    
    expect_equal(results, expected)
})


test_that("CNpreprocessing() must return expected results 03", {
    
    set.seed(2211)
    
    results <- CNpreprocessing(segall=segExample, ratall=rateExample, idCol="ID", 
                               "start", "end", chromCol="chrom", bpStartCol="chrom.pos.start", 
                               bpEndCol="chrom.pos.end",
                               blsize=3, minJoin=0.25, cWeight=0.2, bsTimes=1, chromRange=1, nJobs=1,
                               modelNames="E", normalLength=normalLength,
                               normalMedian=normSegs)
    row.names(results) <- NULL
    
    expected <- segExample
    
    expected$segmedian <- c(0.0866247523, 0.0731923746, 0.0768954424, -0.0073122998, 0.0281455711)
    expected$segmad <- c(0.0640420795, 0.0476423308, 0.0883320200, 0.0410240773, 0.0752964600)
    expected$mediandev <- c(0.071985431920620, 0.058553054211406, 0.062256122026225, -0.021951620158124, 0.013506250751650)
    expected$segerr <- as.double(rep(NA, 5))
    expected$centerz <- c(0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000000, 0.000000000000)
    expected$marginalprob <- c(0.000000000000000, 0.000000000000000, 0.000000000000000, 0.017203615180631, 0.000000000000000)
    expected$maxz <- c(1.000000000000000, 1.000000000000000, 0.538541655965071, 1.0000000000000, 0.990162113892837)
    expected$maxzmean <- c(0.029023619841193, 0.073987616075458, 0.062779292599172, 0.000000000000000, 0.062779292599172)
    expected$maxzsigma <- c(0.002310603657899, 0.002310603657899, 0.002310603657899, 0.011732671852590, 0.002310603657899)
    expected$samplesize <- c(111, 75, 41, 67, 42)
    expected$negtail <- c(1.000000000000000, 1.000000000000000, 1.000000000000000, 0.000372936761926, 0.992841083588778)
    row.names(expected) <- NULL
    
    expect_equal(results, expected)
})

test_that("CNpreprocessing() must return expected results when columns names are not the one by default", {
    
    set.seed(2211)
    
    segExampleTmp <- segExample
    colnames(segExampleTmp) <- c("Id", "Start", "End", "Num.Probes",
                                 "Seg.Median", "Chrom", "Chrom.Pos.Start",
                                 "Chrom.Pos.End", "Cytoband.Start", 
                                 "Cytoband.End", "Abs.Pos.Start", "Abs.Pos.End")
    
    results <- CNpreprocessing(segall=segExampleTmp, ratall=rateExample, idCol="Id", 
                               "Start", "End", chromCol="Chrom", bpStartCol="Chrom.Pos.Start", 
                               bpEndCol="Chrom.Pos.End",
                               blsize=3, minJoin=0.55, cWeight=0.2, bsTimes=1, chromRange=1, nJobs=1,
                               modelNames="E", normalLength=normalLength,
                               normalMedian=normSegs)
    row.names(results) <- NULL
    
    expected <- segExample
    colnames(expected)[1:12] <- c("Id", "Start", "End", "Num.Probes",
                                 "Seg.Median", "Chrom", "Chrom.Pos.Start",
                                 "Chrom.Pos.End", "Cytoband.Start", 
                                 "Cytoband.End", "Abs.Pos.Start", "Abs.Pos.End")
    
    
    expected$segmedian <- c(0.0866247523, 0.0731923746, 0.0768954424, -0.0073122998, 0.0281455711)
    expected$segmad <- c(0.0640420795, 0.0476423308, 0.0883320200, 0.0410240773, 0.0752964600)
    expected$mediandev <- c(0.071985431920620, 0.058553054211406, 0.062256122026225, -0.021951620158124, 0.013506250751650)
    expected$segerr <- as.double(rep(NA, 5))
    expected$centerz <- c(0.0000000000000, 0.0000000000000, 0.0000000000000, 1.0000000000000, 0.0000000000000)
    expected$marginalprob <- c(0.000000000000, 0.000000000000, 0.000000000000, 0.017203615180631, 0.000000000000)
    expected$maxz <- c(1.000000000000, 1.000000000000, 0.538541655965071, 1.000000000000000, 0.990162113892837)
    expected$maxzmean <- c(0.029023619841193, 0.073987616075458, 0.062779292599172, 0.000000000000000, 0.062779292599172)
    expected$maxzsigma <- c(0.002310603657899, 0.002310603657899, 0.002310603657899, 0.011732671852590, 0.002310603657899)
    expected$samplesize <- c(111, 75, 41, 67, 42)
    expected$negtail <- c(1.000000000000, 1.000000000000, 1.000000000000, 0.000372936761926, 0.992841083588778)
    row.names(expected) <- NULL
    
    expect_equal(results, expected)
})




test_that("CNpreprocessing() must return expected results when not ratall", {
    
    results <- CNpreprocessing(segall=segExample, ratall=NULL, "ID", 
                               "start", "end", chromCol="chrom", bpStartCol="chrom.pos.start", 
                               bpEndCol="chrom.pos.end",
                               blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=1, nJobs=1,
                               modelNames="E", normalLength=normalLength,
                               normalMedian=normSegs)
    
    expected <- segExample
    
    expect_equal(results, expected)
})

test_that("CNpreprocessing() must return expected message when not ratall", {
    
    message <- "No raw table, proceeding to comparison"
    
    expect_output(CNpreprocessing(segall=segExample, ratall=NULL, "ID", 
                                "start", "end", chromCol="chrom", bpStartCol="chrom.pos.start", 
                                bpEndCol="chrom.pos.end",
                                blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=1, nJobs=1,
                                modelNames="E", normalLength=normalLength,
                                normalMedian=normSegs), message)
})


test_that("CNpreprocessing() must return error when not idcol", {
    
    message <- "Found unmatched segmented profile IDs"
    
    expect_error(CNpreprocessing(segall=segExample, ratall=rateExample, idCol=NULL, 
                               "start", "end", chromCol="chrom", bpStartCol="chrom.pos.start", 
                               bpEndCol="chrom.pos.end",
                               blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=1,
                               modelNames="E", normalLength=normalLength,
                               normalMedian=normSegs), message)
    
})

test_that("CNpreprocessing() must return expected message when not idcol and no ratall", {
    
    message <- "Found a single segmented profile with no ID \\nNo raw table, proceeding to comparison"
    
    expect_output(CNpreprocessing(segall=segExample, ratall=NULL, idCol=NULL, 
                                 "start", "end", chromCol="chrom", 
                                 bpStartCol="chrom.pos.start", 
                                 bpEndCol="chrom.pos.end",
                                 blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, 
                                 chromRange=1, nJobs=1, modelNames="E", 
                                 normalLength=normalLength,
                                 normalMedian=normSegs), message)
    
})

test_that("CNpreprocessing() must return expected results when not idcol and no ratall", {
    
    results <- CNpreprocessing(segall=segExample, ratall=NULL, idCol=NULL, 
                                    "start", "end", chromCol="chrom", 
                                    bpStartCol="chrom.pos.start", 
                                    bpEndCol="chrom.pos.end",
                                    blsize=5, minJoin=0.25, cWeight=0.4, 
                                    bsTimes=1, chromRange=1, nJobs=1, modelNames="E", 
                                    normalLength=normalLength,
                                    normalMedian=normSegs)
    
    expected <- segExample
    
    expect_equal(results, expected)
})

test_that("CNpreprocessing() must return expected error when not startCol and no bpStartCol", {
    
    message <- "Unable to proceed: incomplete segment annotation"
    
    expect_error(CNpreprocessing(segall=segExample, ratall=rateExample, idCol="ID", 
                               startCol=NULL,  endCol="end", chromCol="chrom", bpStartCol=NULL, 
                               bpEndCol="chrom.pos.end",
                               blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=1,
                               modelNames="E", normalLength=normalLength,
                               normalMedian=normSegs), message)
})

test_that("CNpreprocessing() must return expected error when annotation table given without start and chrom column names", {
    
    message <- paste0("No start and chrom column names provided for annotation table\n")
    
    annotationTmp <- data.frame(PROBEID=c(paste0("ZZ-", 1:10)), CHROM=c(rep(1, 10)), CHROM.POS=c(932544, 1036579, 1212746,
                                                                              1750583, 1750999, 1760583,
                                                                              1850583, 100394039, 101394030,
                                                                              106394039))
    
    expect_error(CNpreprocessing(segall=segExample, ratall=rateExample, idCol="ID", 
                                 startCol=NULL,  endCol=NULL, chromCol="chrom", bpStartCol="chrom.pos.start", 
                                 bpEndCol="chrom.pos.end", annot = annotationTmp,
                                 blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=NULL, nJobs=1,
                                 modelNames="E", normalLength=normalLength,
                                 normalMedian=normSegs), message)
})


test_that("CNpreprocessing() must return expected error when annotation table given without end column name and useEnd set to TRUE", {
    
    message <- paste0("End column name required but not provided in annotation table\n")
    
    annotationTmp <- data.frame(PROBEID=c(paste0("ZZ-", 1:10)), CHROM=c(rep(1, 10)), CHROM.POS=c(932544, 1036579, 1212746,
                                                                                                 1750583, 1750999, 1760583,
                                                                                                 1850583, 100394039, 101394030,
                                                                                                 106394039))
    
    expect_error(CNpreprocessing(segall=segExample, ratall=rateExample, idCol="ID", useEnd=TRUE,
                                 startCol=NULL,  endCol=NULL, chromCol="chrom", bpStartCol="chrom.pos.start", 
                                 bpEndCol="chrom.pos.end", annot = annotationTmp, annotChromCol="CHROM",
                                 annotStartCol="CHROM.POS",annotEndCol=NULL,
                                 blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=NULL, nJobs=1,
                                 modelNames="E", normalLength=normalLength,
                                 normalMedian=normSegs), message)
})

test_that("CNpreprocessing() must return expected error when annotation table given without end column name and useEnd set to FALSE", {
    
    message <- paste0("Incomplete start and end annotation of segments")
    
    annotationTmp <- data.frame(PROBEID=c(paste0("ZZ-", 1:10)), CHROM=c(rep(1, 10)), CHROM.POS=c(932544, 1036579, 1212746,
                                                                                                 1750583, 1750999, 1760583,
                                                                                                 1850583, 100394039, 101394030,
                                                                                                 106394039))
    
    expect_error(CNpreprocessing(segall=segExample, ratall=rateExample, idCol="ID", useEnd=FALSE,
                                    startCol=NULL,  endCol=NULL, chromCol="chrom", bpStartCol="chrom.pos.start", 
                                    bpEndCol="chrom.pos.end", annot = annotationTmp, annotChromCol = "CHROM",
                                    annotStartCol="CHROM.POS",annotEndCol=NULL,
                                    blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=NULL, nJobs=1,
                                    modelNames="E", normalLength=normalLength,
                                    normalMedian=normSegs), message)
})


test_that("CNpreprocessing() must return expected error when not startCol and no annotation table", {
    
    message <- "No annotation table; unable to determine boundary probes/bin"
    
    expect_error(CNpreprocessing(segall=segExample, ratall=rateExample, idCol="ID", 
                                 startCol=NULL,  endCol="end", chromCol="chrom", 
                                 bpStartCol="chrom.pos.start", 
                                 bpEndCol="chrom.pos.end",
                                 blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, 
                                 chromRange=1, nJobs=1, modelNames="E", 
                                 normalLength=normalLength,
                                 normalMedian=normSegs), message)
})


test_that("CNpreprocessing() must return expected error when not idcol and ratall with more than 1 numerical column", {
    
    message <- paste0("Ambiguity: more than 1 numeric column in ",
                      "\"ratall\" matrix. The \"idCol\" ", 
                      "parameter should be specified.\n")
    
    rateExampleTmp <- matrix(data = rep(c(0.0720738404241163, 0.119913919265996, 0.154459489283567, 0.0409946196196852, -0.0828437318817242, 
                                   0.0930527245776551, 0.170908929806745, 0.100289489712059, 0.08662475225992, -0.00385501101587094, -0.195791648504205, 
                                   0.0636341122831496, 0.109449474459845, 0.043428961427158, 0.160174528536138, 0.0410580626539806, 0.0405224616647681, 
                                   0.0986143513758616, 0.0731923745507058, 0.0121261797978406, 0.107795728119586, 0.0832422912347904, 0.0727354592601023, 
                                   0.1475937155592, 0.0810554254709479, 0.162564668451733, 0.0590673628590426, -0.0728672744452223, 0.0284354530744622, 
                                   0.160126815673499, 0.190799021333644, -0.0269711432208357, 0.0203580612494725, -0.00731229981882388, -0.0440263312613861, 
                                   -0.0817496261303279, 0.0146393203392997, 0.0860431917828892, 0.0213604328945665, 0.0789323375946298, 0.0757082196458755, 
                                   -0.0507940506094144, 0.0217577735749916, 0.0281455710909498, 0.0196702655955798, 0.154422985605817, 0.260452283950961, 0.122065762785514, 
                                   -0.0197506768456278, -0.0277870900508885), 2), ncol = 2)
    
    
    expect_error(CNpreprocessing(segall=segExample, ratall=rateExampleTmp, idCol=NULL, 
                                 startCol=NULL,  endCol="end", chromCol="chrom", bpStartCol="chrom.pos.start", 
                                 bpEndCol="chrom.pos.end",
                                 blsize=5, minJoin=0.25, cWeight=0.4, bsTimes=1, chromRange=1, nJobs=1,
                                 modelNames="E", normalLength=normalLength,
                                 normalMedian=normSegs), message)
})
