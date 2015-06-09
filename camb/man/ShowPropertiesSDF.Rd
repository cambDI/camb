\name{ShowPropertiesSDF}
\alias{ShowPropertiesSDF}
\title{Show the fields of a .sdf file
}
\description{
\code{ShowPropertiesSDF} returns a vector containing the names of the fields in an input '.sdf' file.
}
\usage{
ShowPropertiesSDF(structures.file)
}
\arguments{
  \item{structures.file}{
	  A '.sdf' file.
}
}
\value{
A vector containing the field (properties) names present in the input '.sdf' file.
}
\author{
Isidro Cortes-Ciriano <isidrolauscher@gmail.com>
}


\seealso{
\code{\link{GetPropertiesSDF}}
\code{\link{GetPropertySDF}}
}
\examples{
# ShowPropertiesSDF("./test.sdf")
}
