\name{CorrelationPlot}
\alias{CorrelationPlot}
\title{Scatterplot of the Observed against the Predicted Values.
}
\description{
The function depicted a scatterplot of the observed against the predicted values with a machine learning model.
}
\usage{
CorrelationPlot(pred, obs, margin = NULL, main = "", ylab = "Predicted", 
xlab = "Observed", PointSize = 4, ColMargin = "blue", TextSize = 15, 
TitleSize = 15, XAxisSize = 15, YAxisSize = 15, TitleAxesSize = 15, 
tmar = 1, bmar = 1, rmar = 1, lmar = 1, AngleLab = 30, 
LegendPosition = "right", PointColor = "black", PointAlpha = 1, 
PointShape = 16, MarginWidth = 1)
}
\arguments{
  \item{pred}{
  Predicted values.
}
  \item{obs}{
  Observed values.
}
  \item{margin}{
  Bioactivity margin centered in the diagonal of the correlation plot. 
Default value NULL.
}
  \item{main}{
  Plot title.
}
  \item{ylab}{
  Title of the Y axis.
}
  \item{xlab}{
  Title of the X axis.
}
  \item{PointSize}{
  Size of the points.
}
  \item{ColMargin}{
  Color of the bioactivity margin ('margin').
}
  \item{TextSize}{
  Text font size. Default value 15.
}
  \item{TitleSize}{
  Title font size. Default value 15.
}
  \item{XAxisSize}{
  Size of the text on the X axis. Default value 15.
}
  \item{YAxisSize}{
  Size of the text on the Y axis. Default value 15.
}
  \item{TitleAxesSize}{
  Font size of the axes lables. Default value 15.
}
  \item{tmar}{
  Top margin size. Default values is 1.
}
  \item{bmar}{
  Bottom margin size. Default values is 1.
}
  \item{rmar}{
  Right margin size. Default values is 1.
}
  \item{lmar}{
  Left margin size. Default values is 1.
}
  \item{AngleLab}{
  Angle of the labels in the X axis. Default value 30.
}
  \item{LegendPosition}{
  Position of the legend. Default value 'right'.
}
  \item{PointColor}{
  Color of the points in the scatterplot. Default value 'black'.
}
  \item{PointAlpha}{
  Color alpha of the points in the scatterplot. Default value 1.
}
  \item{PointShape}{
  Shape of the points in the scatterplot. Default value 16.
}
  \item{MarginWidth}{
}
}

\value{
A list (ggplot2 plot) with the scatterplot.

}

\author{
Isidro Cortes <isidrolauscher@gmail.com> and Daniel Murrell <dsmurrell@gmail.com>
}

\examples{
CorrelationPlot(pred=seq(1,10)+rnorm(10), obs = seq(1,10),
margin = NULL, main = "", ylab = "Predicted", 
    xlab = "Observed", PointSize = 4, ColMargin = "blue", TextSize = 15, 
    TitleSize = 15, XAxisSize = 15, YAxisSize = 15, TitleAxesSize = 15, 
    tmar = 1, bmar = 1, rmar = 1, lmar = 1, AngleLab = 30, LegendPosition = "right", 
    PointColor = "black", PointAlpha = 1, PointShape = 16, MarginWidth = 1) 
}

