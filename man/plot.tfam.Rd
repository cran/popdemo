% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S3-plot.tfam.R
\name{plot.tfam}
\alias{plot.tfam}
\title{Plot transfer function}
\usage{
\method{plot}{tfam}(x, xvar = NULL, yvar = NULL, mar = c(1.1, 1.1, 0.1, 0.1), ...)
}
\arguments{
\item{x}{an object of class 'tfam' (transfer function analysis matrix) 
created using \code{\link{tfam_lambda}} or \code{\link{tfam_inertia}}.}

\item{xvar, yvar}{(optional) the variables to plot on the x and y axes. May
be \code{"p"}, \code{"lambda"} or \code{"inertia"}. Defaults to
\code{xvar="p"} and \code{yvar="lambda"} for objects created using
\code{tfam_lambda}, and \code{xvar="p"} and \code{yvar="inertia"} for 
objects created using \code{tfam_inertia}.}

\item{mar}{the margin limits on the plots: see \code{\link{par}}}

\item{...}{arguments to be passed to methods: see \code{\link{par}} and
\code{\link{plot}}.}
}
\description{
Plot a matrix of transfer functions
}
\details{
\code{plot.tfam} plots matrices of transfer functions (class 
\code{tfam}) created using \code{\link{tfam_lambda}} or 
\code{\link{tfam_inertia}}. The plot is laid out to correspond with 
the nonzero entries of the matrix used to generate the transfer functions, 
for easy visual comparison of how perturbation affects different matrix 
elements.
}
\examples{
  # Create a 3x3 matrix
  ( A <- matrix(c(0,1,2,0.5,0.1,0,0,0.6,0.6), byrow=TRUE, ncol=3) )

  # Calculate the matrix of transfer functions using default arguments
  ( tfmat<-tfam_lambda(A) )

  # Plot the matrix of transfer functions
  plot(tfmat)

  # Create an initial stage structure
  ( initial <- c(1,3,2) )

  # Calculate the matrix of transfer functions for inertia and 
  # specified initial stage structure using default arguments
  ( tfmat2<-tfam_inertia(A,vector=initial) )

  # Plot the result (defaults to inertia ~ p)
  plot(tfmat2)

  # Plot inertia ~ lambda
  plot(tfmat2, xvar="lambda", yvar="inertia")

}
\seealso{
Constructor functions: \code{\link{tfam_lambda}}, \code{\link{tfam_inertia}}
}
\concept{PVA}
\concept{demography}
\concept{ecology}
\concept{nonlinear}
\concept{perturbation}
\concept{population viability}
\concept{systems control}
\concept{transfer function}
