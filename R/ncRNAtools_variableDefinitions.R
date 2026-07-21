rnaCentralEbiApiURL <- 'https://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral'
rnaCentralApiURL <- 'https://rnacentral.org/api/v1/rna'
rnaCentralRangeSearchURL <- 'https://rnacentral.org/api/v1/overlap/region'

## URLs of the web services used for secondary structure prediction, replacing
## the discontinued ncrna.org/rtools server
viennaRNABaseURL <- 'http://rna.tbi.univie.ac.at/'
ipknotURL <- 'https://ws.sato-lab.org/rtips/ipknot++/'
mcFoldURL <- 'https://www.major.iric.ca/cgi-bin/MC-Fold/mcfold.static.cgi'
rfamBatchSearchURL <- 'https://batch.rfam.org/submit-job'

globalVariables(c("nucleotide1", "nucleotide2", "probability"))

localOS <- Sys.info()["sysname"]

