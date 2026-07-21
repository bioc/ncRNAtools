\name{predictAlternativeSecondaryStructures}
\alias{predictAlternativeSecondaryStructures}
\title{Predicts alternative secondary structures of a given RNA sequence
}
\description{
Identifies potential alternative secondary structures of the provided RNA
sequence with the MC-Fold algorithm, which returns a ranked set of suboptimal
secondary structures. It should be noted that predictions can take considerable
amounts of time, particularly for long sequences.
}
\usage{
predictAlternativeSecondaryStructures(sequence, numStructures=20)
}
\arguments{
  \item{sequence}{string with an RNA sequence whose secondary structure should
  be predicted. Should contain only standard RNA symbols (i.e., "A", "U", "G"
  and "C").
  }
  \item{numStructures}{maximum number of alternative secondary structures to be
  returned, ranked by score. It should be a positive integer. In the default
  behavior, up to 20 structures are returned.
  }
}
\value{
A named list where each element represents a potential secondary structure,
ranked by score, with the first element representing the best-scoring (canonical)
secondary structure. Each element is a list comprising the following two elements:
  \item{sequence }{Query RNA sequence
  }
  \item{secondaryStructure }{Predicted secondary structure
  }
}
\references{
Parisien M, Major F. The MC-Fold and MC-Sym pipeline infers RNA structure from
sequence data. Nature. 2008;452(7183):51-55. doi:10.1038/nature06684

https://www.major.iric.ca/MC-Fold/
}
\examples{
\dontrun{
# Predict alternative secondary structures of an RNA sequence:

alternativeStructures <- predictAlternativeSecondaryStructures("AAAGGGGUUUCCC")

# Count the number of potential alternative structures identified:

length(alternativeStructures)
}
}
