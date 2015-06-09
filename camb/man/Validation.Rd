\name{Validation}
\alias{Validation}
\title{Statistical Metrics to Assess Model Peformance
%%  ~~function to do ... ~~
}
\description{
	The function calculates correlation and error metrics between two 
	numeric vectors.
	These metrics are used to evaluate model performance on a test or external set.
}
\usage{
Validation(pred, obs,resp_tr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pred}{
	  A numeric vector with the predicted values.
}
  \item{obs}{
	  A numeric vector with the observed values.
}
\item{resp_tr}{
	
The argument resp_tr, requires the bioactivity values of the datapoints present in the training set.
These values are required by the metrics Q_{1} and Q_{3} (see camb tutorials).
}
}
\details{
The predictive ability of the models on a test or validation set is evaluated by the calculation of the following statistical metrics:
Q_1^2, 
Q_2^2,
Q_3^2, 
RMSE, 
R^2, and 
R_0^2.
}
\value{
	A list containing the values for the statistical metrics. 
	See the tutorial for the formula of the metrics, 
	and for further details.

}
\references{
camb tutorials.
}
\author{
Isidro Cortes <isidrolauscher@gmail.com> and Daniel Murrell <dsmurrell@gmail.com>
}

\examples{
# See the tutorials for examples.
}
