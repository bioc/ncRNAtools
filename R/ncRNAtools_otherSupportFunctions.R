## Secondary structure prediction backends. The ncrna.org/rtools server that
## previously served CentroidFold, CentroidHomFold, IPknot and RintW was
## discontinued, so each method is now served by a different web service:
##   centroidFold    -> ViennaRNA RNAfold web server
##   centroidHomFold -> Rfam batch search (homology-based, via the matching family)
##   IPknot          -> IPknot web server at the Sato lab
##   RintW           -> MC-Fold (predictAlternativeSecondaryStructures)

sendSecondaryStructureQuery <- function(sequence, method, gammaWeight, inferenceEngine,
                                        alignmentEngine, eValueRfamSearch, numHomSeqsRfamSearch) {
  if (!method %in% c("centroidFold", "centroidHomFold", "IPknot")) {
    stop("Invalid method for secondary structure prediction")
  }
  if (!is.null(inferenceEngine)) {
    checkInferenceEngine(inferenceEngine, method, sequence)
  }
  if (!is.null(gammaWeight)) {
    checkGammaWeight(gammaWeight)
  }
  if (!is.null(alignmentEngine)) {
    checkAlignmentEngine(alignmentEngine, method)
  }
  if (!is.null(eValueRfamSearch)) {
    checkEValueRfamSearch(eValueRfamSearch, method)
  }
  if (!is.null(numHomSeqsRfamSearch)) {
    checkNumHomSeqsRfamSearch(numHomSeqsRfamSearch, method)
  }
  if (method == "centroidFold") {
    return(predictWithViennaRNA(sequence))
  }
  else if (method == "IPknot") {
    return(predictWithIPknot(sequence, inferenceEngine))
  }
  else {
    return(predictWithRfamHomology(sequence))
  }
}

## centroidFold: fold the sequence with the ViennaRNA RNAfold web server, which
## returns a centroid secondary structure and a base pair probability matrix.

predictWithViennaRNA <- function(sequence) {
  message("Running secondary structure prediction. This might take some time.")
  submitURL <- paste(viennaRNABaseURL, "cgi-bin/RNAWebSuite/RNAfold.cgi", sep="")
  response <- POST(submitURL,
                   body=list(PAGE="2", SCREEN=sequence, method="p", proceed=""),
                   encode="multipart")
  responseText <- content(response, as="text", encoding="UTF-8")
  jobID <- sub(".*ID=", "",
               regmatches(responseText, regexpr("PAGE=3&ID=[A-Za-z0-9]+", responseText)))
  if (length(jobID) == 0 || is.na(jobID)) {
    stop("Malformed query or ViennaRNA server unavailable. Please try again.")
  }
  resultsURL <- paste(viennaRNABaseURL, "RNAfold/", jobID, "/", sep="")
  for (attempt in seq_len(120)) {
    Sys.sleep(2)
    if (status_code(GET(paste(resultsURL, "sequence1.vienna", sep=""))) == 200) {
      break
    }
  }
  centroidResponse <- content(GET(paste(resultsURL, "sequence1_centroid.vienna", sep="")),
                              as="text", encoding="UTF-8")
  secondaryStructure <- splitString(splitString(centroidResponse, split="\n")[2], split=" ")[1]
  dotPlotResponse <- content(GET(paste(resultsURL, "sequence1_dp.eps", sep="")),
                             as="text", encoding="UTF-8")
  basePairProbsTable <- parseViennaBasePairProbabilities(dotPlotResponse, sequence)
  return(list(sequence=sequence, secondaryStructure=secondaryStructure,
              basePairProbabilities=basePairProbsTable))
}

## Parse the base pair probabilities from the ViennaRNA dot plot PostScript file,
## where each 'ubox' line holds a base pair as "i j sqrt(probability)", and build
## a table in the format used by generatePairsProbabilityMatrix.

parseViennaBasePairProbabilities <- function(dotPlotResponse, sequence) {
  sequenceCharacters <- splitString(sequence, split="")
  sequenceLength <- length(sequenceCharacters)
  uboxLines <- grep("^[0-9].*ubox$", splitString(dotPlotResponse, split="\n"), value=TRUE)
  pairingProbabilities <- vector(mode="list", length=sequenceLength)
  for (uboxLine in uboxLines) {
    fields <- splitString(trimws(uboxLine), split=" +")
    firstBase <- as.integer(fields[1])
    secondBase <- as.integer(fields[2])
    probability <- as.numeric(fields[3])^2
    pairingProbabilities[[firstBase]] <- c(pairingProbabilities[[firstBase]],
                                           paste(secondBase, format(probability, scientific=FALSE, trim=TRUE), sep=":"))
  }
  maxPairings <- max(c(0, lengths(pairingProbabilities)))
  pairingsMatrix <- matrix("", nrow=sequenceLength, ncol=maxPairings)
  for (position in seq_len(sequenceLength)) {
    if (length(pairingProbabilities[[position]]) > 0) {
      pairingsMatrix[position, seq_along(pairingProbabilities[[position]])] <- pairingProbabilities[[position]]
    }
  }
  basePairProbsTable <- data.frame(Position=seq_len(sequenceLength),
                                   Nucleotide=sequenceCharacters,
                                   pairingsMatrix, stringsAsFactors=FALSE)
  if (maxPairings > 0) {
    colnames(basePairProbsTable)[-(1:2)] <- paste0("Pairing", seq_len(maxPairings))
  }
  return(basePairProbsTable)
}

## IPknot: predict the secondary structure, including pseudoknots, with the
## IPknot web server hosted at the Sato lab.

predictWithIPknot <- function(sequence, inferenceEngine) {
  message("Running secondary structure prediction. This might take some time.")
  serverHandle <- handle(ipknotURL)
  POST(handle=serverHandle, url=paste(ipknotURL, "process.php", sep=""),
       body=list(seq=paste(">", "query", "\n", sequence, sep=""),
                 level="2", model="LinearPartition-C", refinement="1", use_pf="on"),
       encode="multipart")
  viennaResponse <- content(GET(handle=serverHandle, url=paste(ipknotURL, "dload_vienna.php", sep="")),
                            as="text", encoding="UTF-8")
  responseLines <- splitString(viennaResponse, split="\n")
  responseLines <- responseLines[responseLines != ""]
  if (length(responseLines) < 3) {
    stop("Malformed query or IPknot server unavailable. Please try again.")
  }
  return(list(sequence=responseLines[2], secondaryStructure=responseLines[3]))
}

## centroidHomFold: fold the sequence using homology information. The sequence is
## searched against the Rfam database and, if it matches a family, the consensus
## secondary structure of that family is mapped onto the query sequence.

predictWithRfamHomology <- function(sequence) {
  message("Searching the Rfam database for homologs. This might take some time.")
  sequenceFile <- tempfile()
  writeLines(sequence, con=sequenceFile)
  submitResponse <- POST(rfamBatchSearchURL, accept_json(),
                         body=list(sequence_file=upload_file(sequenceFile)),
                         encode="multipart")
  resultURL <- content(submitResponse)$resultURL
  if (is.null(resultURL)) {
    stop("Malformed query or Rfam search server unavailable. Please try again.")
  }
  searchResult <- NULL
  for (attempt in seq_len(100)) {
    Sys.sleep(3)
    searchResult <- content(GET(resultURL, accept_json()))
    if (!is.null(searchResult$closed) && nzchar(searchResult$closed)) {
      break
    }
  }
  if (is.null(searchResult$closed) || !nzchar(searchResult$closed)) {
    stop("The Rfam search did not complete. The service may be temporarily unavailable.")
  }
  if (length(searchResult$hits) == 0) {
    stop("No homologous Rfam family was found for the provided sequence.")
  }
  bestHit <- searchResult$hits[[1]][[1]]
  return(mapRfamConsensusStructure(bestHit, sequence))
}

## Map the consensus secondary structure of the matching Rfam family onto the
## query sequence, keeping the alignment columns that correspond to query
## nucleotides and converting the WUSS annotation to Dot-Bracket notation.

mapRfamConsensusStructure <- function(hit, sequence) {
  alignedStructure <- splitString(hit$alignment$ss, split=" +")[2]
  alignedQuery <- splitString(hit$alignment$user_seq, split=" +")[3]
  structureCharacters <- splitString(alignedStructure, split="")
  queryCharacters <- splitString(alignedQuery, split="")
  queryColumns <- grepl("[AUGCaugc]", queryCharacters)
  wussStructure <- structureCharacters[queryColumns]
  dotBracketStructure <- ifelse(wussStructure %in% c("(", "[", "{", "<"), "(",
                                ifelse(wussStructure %in% c(")", "]", "}", ">"), ")", "."))
  return(list(sequence=sequence,
              secondaryStructure=paste(dotBracketStructure, collapse="")))
}

## RintW: predict a set of alternative secondary structures with MC-Fold, which
## returns a ranked set of suboptimal structures for the provided sequence.

sendAlternativeSecondaryStructureQuery <- function(sequence, gammaWeight, inferenceEngine, numStructures) {
  message("Running secondary structure prediction. This might take some time.")
  response <- GET(mcFoldURL,
                  query=list(pass="lucy", sequence=sequence, top=as.character(numStructures),
                             explore="15", name="", mask=""),
                  timeout(300))
  responseText <- content(response, as="text", encoding="latin1")
  responseText <- gsub("<[^>]*>", "", responseText)
  structureLines <- grep("^[0-9]+\\)[[:space:]]*[().]+",
                         trimws(splitString(responseText, split="\n")), value=TRUE)
  alternativeStructures <- sub("^[0-9]+\\)[[:space:]]*([().]+).*", "\\1", trimws(structureLines))
  alternativeStructures <- unique(alternativeStructures)
  if (length(alternativeStructures) == 0) {
    stop("Malformed query or MC-Fold server unavailable. Please try again.")
  }
  return(list(sequence=sequence, alternativeStructures=alternativeStructures))
}

retrieveAlternativeSecondaryStructureResults <- function(queryResult) {
  alternativeStructures <- lapply(queryResult$alternativeStructures,
                                  function(structure) list(sequence=queryResult$sequence,
                                                           secondaryStructure=structure))
  names(alternativeStructures) <- c("alternativeStructure1-Canonical",
                                    if (length(alternativeStructures) > 1) {
                                      paste("alternativeStructure", seq(2, length(alternativeStructures)), sep="")
                                    })
  return(alternativeStructures)
}

extractBasePairProbability <- function(basePairProbsTableField) {
  if (basePairProbsTableField != "") {
    basePairProb <- as.numeric(splitString(basePairProbsTableField, split=":"))
    names(basePairProb) <- c("targetPosition", "probability")
    return(basePairProb)
  }
}

makeCompositeMatrix <-function(basePairProbsMatrix, pairedBases) {
  compositeMatrix <- basePairProbsMatrix
  compositeMatrix[upper.tri(compositeMatrix)] <- 0
  for (row in seq_len(nrow(pairedBases))) {
    compositeMatrix[pairedBases[row, "Position1"], pairedBases[row, "Position2"]] <- 1
  }
  return(compositeMatrix)
}

flattenDotBracket <- function(extendedDotBracketString) {
  extendedDotBracketVector <- splitString(extendedDotBracketString, "")
  simpleDotBracketVector <- character(length=length(extendedDotBracketVector))
  for (i in seq_len(length(extendedDotBracketVector))) {
    if (extendedDotBracketVector[i] %in% c("(", "[", "<", "{", "A", "B", "C", "D")) {
      simpleDotBracketVector[i] <- "("
    } else if (extendedDotBracketVector[i] %in% c(")", "]", ">", "}", "a", "b", "c", "d")) {
      simpleDotBracketVector[i] <- ")"
    } else if (extendedDotBracketVector[i] == ".") {
      simpleDotBracketVector[i] <- "."
    } else {
      stop("Invalid characters present in extended Dot-Bracket string")
    }
  }
  return(paste(simpleDotBracketVector, collapse=""))
}

splitString <- function(string, split, ...) {
  return(unlist(strsplit(string, split=split, ...)))
}

removeNewLines <- function(string) {
  return(gsub("[\r\n]", "", string))
}
