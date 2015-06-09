\name{MaxPerf}
\alias{MaxPerf}
\title{Distribution of Maximum Theoretical Values Achievable given the Dataset and its Uncertainty
}
\description{
	Calculates the ditribution of model validation metrics that are achievable given the size of the dataset,
	the uncertainty in the response variable, and the distribution of the responsible variable quantified by its mean and standard deviation.
	Therefore, these distributions help to assess models overfitting; e.g. a model trained on a dataset with high uncertainty exhibiting
	high correlation values might be overoptimistic.
	See the tutorials for more information and examples.
}
\usage{
MaxPerf(meanNoise = 0, sdNoise, resp, lenPred, stds = NULL,     iters = 1000, filename = NULL, pdfW = 10, pdfH = 10, TextSize = 15,     TitleSize = 15, XAxisSize = 15, YAxisSize = 15, TitleAxesSize = 15,     tmar = 1, bmar = 1, rmar = 1, lmar = 1, AngleLab = 30, 
    LegendPosition = "right")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{meanNoise}{
Mean value of the noise in the data. Default value 0.
}
  \item{sdNoise}{
Standard deviation of the noise in the data. See the work by Kramer et al. about uncertainty in public bioactivity databases.
}
  \item{resp}{
Vector containing the values for the dependent variable in the dataset (y).
}
  \item{lenPred}{
Number of datapoints of the external (hold-out) set.
%%     ~~Describe \code{lenPred} here~~
}
\item{stds}{
Vector containing the experimental error for the response variable. The default value is NULL.
}
  \item{iters}{
	  Number of iterations. Default value 1000.
%%     ~~Describe \code{iters} here~~
}
  \item{filename}{
If not NULL, file where the plot will be saved. Default value NULL.
%%     ~~Describe \code{filename} here~~
}
  \item{pdfW}{
Width of the .pdf file, in centimeters, where the plot will be saved. Default value 10.
%%     ~~Describe \code{pdfW} here~~
}
  \item{pdfH}{
Height of the .pdf file, in centimeters, where the plot will be saved. Default value 10.
%%     ~~Describe \code{pdfH} here~~
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
%%     ~~Describe \code{XAxisSize} here~~
Fontsize of the X axis. Default value 15.
}
  \item{YAxisSize}{
%%     ~~Describe \code{YAxisSize} here~~
Fontsize of the Y axis. Default value 15.
}
  \item{TitleAxesSize}{
%%     ~~Describe \code{TitleAxesSize} here~~
Fontsize of the axes titles. Default value 15.
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
	  Angle of the labels of the X axis. Default value 30.
%%     ~~Describe \code{AngleLab} here~~
}
  \item{LegendPosition}{
	  Position of the legend. Default value 'right'.
%%     ~~Describe \code{LegendPosition} here~~
}
}
\references{
Cortes-Ciriano et al. 'Proteochemometric Modeling in a Bayesian Framework'. J. Cheminf. 6, 35. 2014
\url{http://www.jcheminf.com/content/6/1/35}
}
\author{
Isidro Cortes <isidrolauscher@gmail.com> and Daniel Murrell <dsmurrell@gmail.com>
}



