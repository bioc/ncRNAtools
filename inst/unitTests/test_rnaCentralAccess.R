library(RUnit)
library(GenomicRanges)
library(IRanges)
library(S4Vectors)

## The tests below depend on the RNAcentral web service and are only run when it
## returns valid results, so that transient outages do not cause spurious errors.

## Test rnaCentralTextSearch
textSearchResult <- tryCatch(rnaCentralTextSearch("HOTAIR"), error=function(e) character(0))
if (length(textSearchResult) > 0) {
  checkTrue(grepl("^URS", textSearchResult[1]))
}

## Test rnaCentralRetrieveEntry
retrievedEntry <- tryCatch(rnaCentralRetrieveEntry("URS000075C808_9606"), error=function(e) NULL)
if (!is.null(retrievedEntry) && !is.null(retrievedEntry$rnaCentralID)) {
  checkEquals("URS000075C808_9606", retrievedEntry$rnaCentralID)
}

## Test rnaCentralGenomicCoordinatesSearch
genomicSearchResult <- tryCatch(
  rnaCentralGenomicCoordinatesSearch(GenomicRanges::GRanges(seqnames=S4Vectors::Rle("chr3"),
                                                            ranges=IRanges::IRanges(39788104, 39795326)),
                                     "Homo sapiens"),
  error=function(e) NULL)
if (!is.null(genomicSearchResult) && length(genomicSearchResult) > 0) {
  checkTrue(grepl("^URS", genomicSearchResult[[1]][[1]]$rnaCentralID))
}
