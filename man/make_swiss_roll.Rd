% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tsne_fit.r
\name{make_swiss_roll}
\alias{make_swiss_roll}
\title{Swiss roll data set}
\usage{
make_swiss_roll(nn)
}
\arguments{
\item{nn}{The number of data points to be generated.}
}
\value{
A data frame with four variables. The three first are the swiss
roll data points, the fourth is vector fo colors useful for plotting.
}
\description{
Simulate swiss roll data. The swiss roll data set is used to test various dimension reduction algorithms.
}
\section{Reference}{

Dinoj Surendran, Swiss Roll Dataset
\url{http://people.cs.uchicago.edu/~dinoj/manifold/swissroll.html}
}

\examples{
xx <- make_swiss_roll(150)
pairs(xx[,1:3], col=xx[,4], pch=16)

}
