## Local Whittle Estimation Method
##------------------------------------------------------------------------
## Original Code written by: Fotis Papailias
## based on and refining parts of the original code written by Taqqu, M. S.
## for S-PLUS.
##------------------------------------------------------------------------
## Notes: This method introduced by Robinson (1995) it is very convenient
## providing good estimation results with relatively small MSE and bias.
## However, this estimator can *not* be used for locally stationary or
## non-stationary series; is must be used STRICTLY on
## stationary Long Memory Time Series.
##------------------------------------------------------------------------
## Requirements: (longmemo) package
##------------------------------------------------------------------------
## References:
##
## [1]
##  Robinson, P. M. (1995): “Gaussian Semi-parametric Estimation of
##                            Long Range Dependence”,
##  Annals of Statistics, 23, 1630–1661.
##
## [2]
## http://math.bu.edu/people/murad/methods/locwhitt/locwhitt.ps
##------------------------------------------------------------------------
## Notes: "d" is the long memory parameter, "x" is the series
## and check [2] of the above references to find out the "im"
## about how to define the lag. "im" is the parameter that
## defines the lag as m = length of series/ im
##------------------------------------------------------------------------
## You can call it using "lwhittle(x, im)" where "x" the series and
## "im" the parameter we described above.
##------------------------------------------------------------------------

require("longmemo") # only for per()  [[ MM: FIXME! ]]
lw <- function(d, x, im)
{
    h <- d+0.5
    k <- length(peri1 <- per(x))
    len <- length(x)
    m <- floor(len / im)
    if(m+1 > k)
        stop("'im' must be larger than ", format(len/(k-1)))
    peri <- peri1[2:(m+1)]
    z <- 1:m
    freq <- ((2*pi)/len) * z
    ## result
    log(sum(freq^(2*h-1)*peri)) - (2*h)/m * sum(log(freq))
}

lwhittle <- function(x, im)
{
    optimize(lw, interval = c(-.5, 1), ## d = h - 1/2
             x, im)
}
