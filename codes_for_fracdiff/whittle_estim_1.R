## Whittle's Likelihood Estimation Method #1
##------------------------------------------------------------------------
## Original Code written by: Wilfredo Palma
## http://www.mat.puc.cl/~wilfredo/english/html/theory.html
##------------------------------------------------------------------------
## Edited by: Fotis Papailias (f.papailias@qmul.ac.uk)
## Notes: This method works fine for the 'nice' ARFIMA (0, d, 0)
## 	model and in addition provides (some) good results for other
##	models too. More demanding users should use Method #2 instead.
##------------------------------------------------------------------------
## Reference: (Lots. A useful review is:)
##	Taqqu, M. S., Teverovsky, V. and W. Willinger (1995):
##	“Estimators for Long-Range Dependence: An Empirical Study”,
##	Fractals, 3 (4), 785-802.
##------------------------------------------------------------------------
## http://math.bu.edu/people/murad/pub/estimators-posted.ps
##------------------------------------------------------------------------
## You can call the function, in general, as "whittle1(x)" where "x" is the
## underlying time series.
##------------------------------------------------------------------------

###  MM: FIXME -- very similar code
###  has been hidden in package  ../../longmemo/R/WhittleEst.R !!

whittle1 <- function(series)
{
### MM: FIXME! -- use optimize for 1-dim optimization

    optim(par = 0.25, whittle.loglik, gr = NULL,
          method = "L-BFGS-B", series = series)$par

}

##------------------------------------------------------------------------
## AUXILIARY FUNCTIONS
##------------------------------------------------------------------------

whittle.loglik <- function(x, series)
{
    fn.density <- function(x, d)
    {
        ((2 * sin(x/2))^(-2 * d)) / (2*pi)
    }

    series <- series - mean(series)
    a <- fft(series)
    a <- Mod(a)^2
    n <- length(series)
    a <- a/(2 * pi * n)
    m <- n/2
    w <- (2 * pi * (1:m))/n
    b <- fn.density(w, x)
    sigma2 <- (2 * sum(a[1:m]/b))/n
    loglik <- 2 * pi * (sum(log(b)) + sum(a[1:m]/b)/sigma2)
    return(loglik/n + pi * log(sigma2))
}

