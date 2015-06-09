\name{PairwiseDist}
\alias{PairwiseDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pairwise Distance (Similarity) Matrix
%%  ~~function to do ... ~~
}
\description{
	The function is based on the vegdist function from the vegan package.
	It calculated the pairwise distance similarity matrix for all vectors
	input in a matrix or data.frame.
	The functions operates on a row basis.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PairwiseDist(Data, method = "jaccard", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
	  A numeric data.frame or matrix containing compound, 
	  protein or amino acid descriptors (or any combination thereof).
%%     ~~Describe \code{Data} here~~
}
  \item{method}{
Available distance metrics are: "manhattan", "euclidean", "canberra", "bray", 
	            "kulczynski", "jaccard", "gower", "altGower", "morisita", 
				            "horn", "mountford", "raup", "binomial", "chao", 
							            "cao".
See the documentation of the R package vegan for details.
}
}
\details{
	For further detials see the documentation in the R package vegan.
}
\value{
	A data.frame with the all pairwise distances.
}

\author{
Isidro Cortes <isidrolauscher@gmail.com> and Daniel Murrell <dsmurrell@gmail.com>
}

\seealso{
	PairwiseDistPlot
}
\examples{
	m = matrix(abs(rnorm(20)),4,4) 
	mDist = PairwiseDist(m)
	head(mDist)
}

