rnaCentralEbiApiURL <- 'https://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral'
rnaCentralApiURL <- 'https://rnacentral.org/api/v1/rna'
rnaCentralRangeSearchURL <- 'https://rnacentral.org/api/v1/overlap/region'

rtoolsBaseURL <- 'https://www.ncrna.org/rtools/'

globalVariables(c("nucleotide1", "nucleotide2", "probability"))

localOS <- Sys.info()["sysname"]

