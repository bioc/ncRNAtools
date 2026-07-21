library(RUnit)

testSequence <- "GGGGCUAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCAUAGCUCCACCA"
testStructure <- "((((((((.(((((((((..(((.(.(((((.......))))).).))).))))).))....))))))))))...."
testBasePairProbabilitiesTable <- read.csv(system.file("extdata", "exampleBasePairProbabilitiesTable.csv", package="ncRNAtools"))
testPairedBases <- findPairedBases(testStructure, testSequence)

## Test findPairedBases

checkTrue(is.data.frame(testPairedBases))

## Test pairsToSecondaryStructure

checkEquals(pairsToSecondaryStructure(testPairedBases, testSequence), testStructure)

## Test generatePairsProbabilityMatrix

checkTrue(is.matrix(generatePairsProbabilityMatrix(testBasePairProbabilitiesTable)))

## Tests of the prediction functions, which depend on external web services and
## are only run when those services are reachable

ipknotReachable <- tryCatch({
  httr::HEAD("https://ws.sato-lab.org/rtips/ipknot++/", httr::timeout(15))
  TRUE
}, error=function(e) FALSE)

if (ipknotReachable) {
  ## Test predictSecondaryStructure
  ipknotPrediction <- predictSecondaryStructure(testSequence, "IPknot")
  checkEquals(nchar(ipknotPrediction$secondaryStructure), nchar(testSequence))
}
