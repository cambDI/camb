\name{PairwiseDistPlot}
\alias{PairwiseDistPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distribution of Pairwise Similarities
%%  ~~function to do ... ~~
}
\description{
	The function depicts the distribution of pairwise similarities.
}
\usage{
PairwiseDistPlot(Data, xlab = "", ylab = "", main = "", TextSize = 15, 
TitleSize = 15, XAxisSize = 15, YAxisSize = 15, TitleAxesSize = 15, tmar = 1, 
bmar = 1, rmar = 1, lmar = 1, AngleLab = 30, binwidth = NULL, fillCol = "white",
 Colour = "black", DensityFill = "#FF6666", DensityAlpha = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
A data.frame with a single column named 'Distance'. This is the default output of PairwiseDist.

%%     ~~Describe \code{Data} here~~
}
  \item{xlab}{
	  Label of the X axis.
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
	  Label of the Y axis.
%%     ~~Describe \code{ylab} here~~
}
  \item{main}{
Plot title.
%%     ~~Describe \code{main} here~~
}
  \item{TextSize}{
	  Fontsize of the text in the plot. Default value 15.
%%     ~~Describe \code{TextSize} here~~
}
  \item{TitleSize}{
	  Fontsize of the title. Default value 15.
%%     ~~Describe \code{TitleSize} here~~
}
  \item{XAxisSize}{

	  Fontsize of the X axis. Default value 15.
%%     ~~Describe \code{XAxisSize} here~~
}
  \item{YAxisSize}{
	  Fontsize of the Y axis. Default value 15.
%%     ~~Describe \code{YAxisSize} here~~
}
  \item{TitleAxesSize}{
	  Fontsize of both the X and Y axes. Default value 15.
%%     ~~Describe \code{TitleAxesSize} here~~
}
  \item{tmar}{
	  Top margin size. Default value 1.
%%     ~~Describe \code{tmar} here~~
}
  \item{bmar}{
	  Bottom margin size. Default value 1.
%%     ~~Describe \code{bmar} here~~
}
  \item{rmar}{
	  Right margin size. Default value 1.
%%     ~~Describe \code{rmar} here~~
}
  \item{lmar}{

	  Left margin size. Default value 1.
%%     ~~Describe \code{lmar} here~~
}
  \item{AngleLab}{
	  Angle of the labels in the X axis. Default value 30.
%%     ~~Describe \code{AngleLab} here~~
}
  \item{binwidth}{
	  Width of the bins of the hitogram. Default NULL, 
	  which corresponds to 1/30 of the range of the data (see ??stat_bin of ggplot2).
%%     ~~Describe \code{binwidth} here~~
}

  \item{fillCol}{
	  Fill color of the histogram. Default value 'white'.
%%     ~~Describe \code{fillCol} here~~
}
  \item{Colour}{
	  Line color of the histogram. Default value 'black'.
%%     ~~Describe \code{Colour} here~~
}
  \item{DensityFill}{
	  Fill color of the distribution. Default value "#FF6666".
%%     ~~Describe \code{DensityFill} here~~
}
  \item{DensityAlpha}{
	  Alpha for the fill color of the distribution. Default value 0.2.
}
}
\value{
	A ggplot2 object with the pairwise distance distribution.
}
\references{
Package ggplot2.
}
\author{
Isidro Cortes <isidrolauscher@gmail.com> and Daniel Murrell <dsmurrell@gmail.com>
}



\seealso{
	PairwiseDist
}
\examples{
m <- matrix(abs(rnorm(1600)),40,40)
mDist <- PairwiseDist(m)
head(mDist)
mDistPlot <- PairwiseDistPlot(mDist,xlab = "", ylab = "", main = "", TextSize = 15, 
TitleSize = 15, XAxisSize = 15, YAxisSize = 15, TitleAxesSize = 15, tmar = 1, bmar = 1,
 rmar = 1, lmar = 1, AngleLab = 30, binwidth = NULL, fillCol = "white",
  Colour = "black", DensityFill = "#FF6666", DensityAlpha = 0.2)


}

