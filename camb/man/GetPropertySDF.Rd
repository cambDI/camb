\name{GetPropertySDF}
\alias{GetPropertySDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract the values corresponding to a given property for the molecules from an input '.sdf' file.
}
\description{
	This function returns a data.frame containing the values for a selected property for
	a user-defined number of molecules from an input '.sdf' file.

}
\usage{
GetPropertySDF(structures.file, property = "", number_processed = -1)
}
\arguments{
  \item{structures.file}{
	  Input '.sdf' file.
}
  \item{property}{
Property name (in string format) which values should be returned.
}
  \item{number_processed}{
	  The number of molecules for which the propertt values should be returned.
	  The property values for all molecules are returned when using the default value, namely -1.

}
}
\value{
A data.frame containing the value for the selected property (argument \code{property})
}
\author{
Isidro Cortes-Ciriano <isidrolauscher@gmail.com>}

\seealso{
\code{\link{GetPropertiesSDF}}
\code{\link{ShowPropertiesSDF}}
}
