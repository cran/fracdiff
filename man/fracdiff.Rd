\name{fracdiff}
\title{fracdiff: Maximum likelihood parameter estimates for}
          fractionally-differenced ARIMA (p,d,q) models

\usage{
fracdiff( x, nar = 0, nma = 0, dtol = <see below>, M = 100)

}
\arguments{
\item{x}{
  time series for the ARIMA model}
\item{nar}{
  number of autoregressive parameters}
\item{nma}{
  number of moving average parameters}
\item{dtol}{
  interval of uncertainty for d
  If dtol is less than zero,
  the fourth root of machine precision will be used.
  dtol will be altered if necessary by the program.}
\item{M}{
  number of terms in the likelihood approximation (see Haslett and
  Raftery 1989)}
}
\description{
Calculates the maximum likelihood estimators of the parameters
of a fractionally-differenced ARIMA (p,d,q) model, together (if possible)
with their estimated covariance and correlation matrices and
standard errors, as well as the value of the maximized likelihood.
The likelihood is approximated using the fast and accurate method
of Haslett and Raftery (1989).
}
\value{
a list containing the following elements :
\item{log.likelihood}{
  logarithm of the maximum likelihood}
\item{d}{
  optimal fractional-differencing parameter}
\item{ar}{
  vector of optimal autoregressive parameters }
\item{ma}{
  vector of optimal moving average parameters}
\item{covariance.dpq}{
  covarianvce matrix of the parameter estimates 
  (order : d, ar, ma)}
\item{stderror.dpq}{
  standard errors of the parameter estimates 
  (order : d, ar, ma)}
\item{correlation.dpq}{
  correlation matrix of the parameter estimates 
  (order : d, ar, ma)}
\item{dtol}{
  interval of uncertainty for d}
}
\note{
Ordinarily nar and nma should not be too large (say < 10) 
to avoid degeneracy in the model. 
The function \code{fracdiff.sim} is available for generating test problems.
}
\section{Method}{
The optimization is carried out in two levels : an outer univariate unimodal
optimization in d over the interval [0,.5] (uses Brent's fmin algorithm), and
an inner nonlinear least-squares optimization in the AR and MA parameters to
minimize white noise variance (uses the MINPACK subroutine \code{lm}DER).
written by Chris Fraley (March 1991)
}
\references{
J. Haslett and A. E. Raftery, "Space-time Modelling with Long-memory
Dependence: Assessing Ireland's Wind Power Resource (with Discussion)",
Applied Statistics, 38, 1-50.

R. Brent, Algorithms for Minimization without Derivatives, Prentice-Hall
(1973). 

J. J. More, B. S. Garbow, and K. E. Hillstrom, 
Users Guide for MINPACK-1, Technical Report ANL-80-74, 
Applied Mathematics Division, Argonne National Laboratory (August 1980).
}
\seealso{
  \code{fracdiff.sim}
}
\examples{
ts.test <- fracdiff.sim( 5000, ar = .2, ma = -.4, d = .3)
fracdiff( ts.test$series, nar = length(ts.test$ar), nma = length(ts.test$ma))
}



