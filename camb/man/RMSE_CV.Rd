\name{RMSE_CV}
\alias{RMSE_CV}
\title{Extract the cross validated RMSE (or another metric) from a caret model}
\usage{
RMSE_CV(model, digits = 3, metric = "RMSE")
}
\arguments{
	\item{mode}{
		Trained caret model from which the value for the cross-validated metric specified in metric.
	}
	\item{digits}{
	Number of decimals to which the value is to be rounded.
	}
	\item{metric}{
	Cross-validated metric which value is to be obtained. Default value is 'RMSE'.
	}
	
	}


\description{
Extract the value of a cross-validated metric (e.g. RMSE) from a caret model
}
\author{
Isidro Cortes <isidrolauscher@gmail.com> and Daniel Murrell <dsmurrell@gmail.com>
}
