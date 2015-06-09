\name{AADescs}
\alias{AADescs}
\title{Amino Acid Descriptor Calculation}
\usage{
AADescs(Data, type = "Z5", ...)
}
\arguments{
\item{Data}{A character, vector, matrix or data.frame containing the amino acids in either one-letter or three-letter format.
Amino acids symbols are valid in capitals or in lower-case.
Gaps in alignments are expected to be represented with the charachter "-".
The value of any of the descriptors provided for "-" is 0.}

\item{type}{Type of descriptors to be calculated. Default value is 5 z-scales.
Any combination of descriptors is valid.
A vector containing the abbreviation of the desired descriptors is taken as argument.}
}
\value{
A data.frame which columns are indexed by the descriptors, and
rows by the rows amino acid sequence.
If the input data is a matrix or data.frame, the number of rows in the original data.frame or matrix,
and the number of rows of the ouput data.frame are equal.
If several descriptor types are chosen, descriptors are concatenated for the ease of further modeling.
Column names indicate the amino acid position in the original input data, and the type of descriptor.
}
\description{
The function calculates amino acid descriptors for natural amino acids.
Currently available descriptors are: 3 and 5 z-scales ("Z3" and "Z5"),
T-scales ("TScales"), ST-scales ("STScales"),
principal components score Vectors of Hydrophobic, Steric, and Electronic properties ("VHSE"),
BLOSUM ("BLOSUM"), FASGAI ("FASGAI"), MSWHIM ("MSWHIM"), and ProtFP ("ProtFP8").
See references for further information on these descriptors.
}
\examples{
AADescs(c("A","A"))
}
\author{
Isidro Cortes <isidrolauscher@gmail.com> and Daniel Murrell <dsmurrell@gmail.com>
}
\references{
\url{http://www.jcheminf.com/content/5/1/41}

\url{http://www.jcheminf.com/content/5/1/42}
}

