\name{GetPropertiesSDF}
\alias{GetPropertiesSDF}
\title{Extract the properties values for the molecules in an input '.sdf' file.
}
\description{
This function returns a data.frame containing the values for all properties for
a user-defined number of molecules from an input '.sdf' file.
}
\usage{
GetPropertiesSDF(structures.file, number_processed = -1)
}
\arguments{
  \item{structures.file}{
Input '.sdf' file.
}
  \item{number_processed}{
	  The number of molecules which properties should be returned.
	  The properties values for all molecules are returned when using the default value, namely -1.
}
}
\value{
A data.frame containing the values for the properties.
}
\author{
Isidro Cortes-Ciriano <isidrolauscher@gmail.com>
}


\seealso{
\code{\link{ShowPropertiesSDF}}
\code{\link{GetPropertySDF}}
}
