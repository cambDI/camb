\name{ErrorBarplot}
\alias{ErrorBarplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Barplot with error bars.
%%  ~~function to do ... ~~
}
\description{
	The function creates easily customizable barplots with error bars.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ErrorBarplot(X, Y, err, fill = X, main = "", ylab = "", xlab = "", TextSize = 15, TitleSize = 15, XAxisSize = 15, YAxisSize = 15, TitleAxesSize = 15, AngleLab = 35, barcol = "red", barSize = 1, barWidth = 0.3, LegendName = "Legend", ColLegend = 1, RowLegend = NULL, LegendPosition = "right", tmar = 1, bmar = 1, rmar = 1, lmar = 1, stat = "identity")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	  A vector containing the aesthetics corresponding to the X axis.

%%     ~~Describe \code{X} here~~
}
  \item{Y}{
	  The values for the ordenate axis.
%%     ~~Describe \code{Y} here~~
}
  \item{err}{
	  The standard deviation corresponding to the Y values.
%%     ~~Describe \code{err} here~~
}
  \item{fill}{
	  The groups that will be used to color the bars.
	  The groups defined in the X variable are used by default.
%%     ~~Describe \code{fill} here~~
}
  \item{main}{
	  Plot title.
%%     ~~Describe \code{main} here~~
}
  \item{ylab}{
	   Title of the Y axis.
%%     ~~Describe \code{ylab} here~~
}
  \item{xlab}{
	   Title of the X axis.
	   
%%     ~~Describe \code{xlab} here~~
}
  \item{TextSize}{
	  Text font size. Default value 15.

%%     ~~Describe \code{TextSize} here~~
}
  \item{TitleSize}{
%%     ~~Describe \code{TitleSize} here~~
Title font size. Default value 15.
}
  \item{XAxisSize}{
	   Size of the text on the X axis. Default value 15.
%%     ~~Describe \code{XAxisSize} here~~
}
  \item{YAxisSize}{
	  Size of the text on the Y axis. Default value 15.
%%     ~~Describe \code{YAxisSize} here~~
}
  \item{TitleAxesSize}{
	  Size of the title of both the X and Y axis. Default value 15.
%%     ~~Describe \code{TitleAxesSize} here~~
}
  \item{AngleLab}{
%%     ~~Describe \code{AngleLab} here~~
Angle of the labels in the X axis. Default value 30.
}
  \item{barcol}{
	  Color of the error bars. Default value "red".
%%     ~~Describe \code{barcol} here~~
}
  \item{barSize}{
	  Size of the error bars. Default value 1.
%%     ~~Describe \code{barSize} here~~
}
  \item{barWidth}{
	  Width of the error bars. Default value 0.3.
%%     ~~Describe \code{barWidth} here~~
}
  \item{LegendName}{
%%     ~~Describe \code{LegendName} here~~
Name of the legend.
}
  \item{ColLegend}{
	  Number of columns of the legend. Default value 1.
%%     ~~Describe \code{ColLegend} here~~
}
  \item{RowLegend}{
	   Number of rows of the legend. Default value NULL - as many as groups-.
%%     ~~Describe \code{RowLegend} here~~
}
  \item{LegendPosition}{
	  Position of the legend. Default value "right".
%%     ~~Describe \code{LegendPosition} here~~
}
  \item{tmar}{
	  Top margin size. Default values is 1.

%%     ~~Describe \code{tmar} here~~
}
  \item{bmar}{
	  Bottom margin size. Default values is 1.
%%    ~~Describe \code{bmar} here~~
}
  \item{rmar}{
%%     ~~Describe \code{rmar} here~~
Right margin size. Default values is 1.
}
  \item{lmar}{
	   Left margin size. Default values is 1.
%%     ~~Describe \code{lmar} here~~
}
  \item{stat}{
	  Default value "identity".
}
}
\value{
	A list containing the barplot (ggplot2 object).
}
\author{
Isidro Cortes <isidrolauscher@gmail.com> and Daniel Murrell <dsmurrell@gmail.com>
}


\examples{
	d = data.frame(Y=seq(1,4),err=rep(1,4),X=c("A","B","C","D"))

# Example 1
ErrorBarplot(d$X,d$Y,d$err,fill=d$X,
main = "", ylab = "", xlab = "", 
TextSize = 15, TitleSize = 15, XAxisSize = 15, YAxisSize = 15, 
TitleAxesSize = 15, AngleLab = 35, barcol = "red", barSize = 1, 
barWidth = 0.3, LegendName = "Legend", ColLegend = 1, RowLegend = NULL, 
LegendPosition = "right", tmar = 1, bmar = 1, rmar = 1, lmar = 1, 
stat = "identity")


# Example 2
ErrorBarplot(d$X,d$Y,d$err,fill=d$X,
main = "Example 2 ErrorBarplot", ylab = "Value", xlab = "Group",
TextSize = 15, TitleSize = 15, XAxisSize = 15, YAxisSize = 15, 
TitleAxesSize = 15, AngleLab = 0, barcol = "green", barSize = 1, 
barWidth = 0.6, LegendName = "Example Legend", ColLegend = 1, 
RowLegend = NULL, LegendPosition = "right", 
tmar = 1, bmar = 1, rmar = 1, lmar = 1, 
stat = "identity")


}

