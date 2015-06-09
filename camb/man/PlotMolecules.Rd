\name{PlotMolecules}
\alias{PlotMolecules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Compounds from a .sdf File.
%%  ~~function to do ... ~~
}
\description{
	The function plots the chemical structures provided in a .sdf file.
	The plots can also be sent to a 2x2 grid in a .pdf file.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PlotMolecules(sdf.file, IDs, pdf.file = NULL, PDFMain = NULL, useNameAsTitle = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sdf.file}{
	  The .sdf file with the molecules.
%%     ~~Describe \code{sdf.file} here~~
}
  \item{IDs}{
	  The IDs of the molecules to be depicted (the ordinal position of the molecules in the .sdf file). 
	  Currently, a maximum of four IDs is supported.
%%     ~~Describe \code{IDs} here~~
}
  \item{pdf.file}{
	  If not NULL, the .pdf where the molecules will be depicted.
%%     ~~Describe \code{pdf.file} here~~
}
  \item{PDFMain}{
If not NULL, the title of the molecule depiction in the .pdf file.
%%     ~~Describe \code{PDFMain} here~~
}
  \item{useNameAsTitle}{
	  If TRUE, the names of the molecules as especified in the .sdf file are used as molecules names in the depiction.
}
}
\value{
	A list with the plots of the molecules.
}
\author{
Isidro Cortes <isidrolauscher@gmail.com> and Daniel Murrell <dsmurrell@gmail.com>
}




