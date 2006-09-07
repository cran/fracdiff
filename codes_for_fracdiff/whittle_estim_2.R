## Whittle's Likelihood Estimation Method #2
##------------------------------------------------------------------------
## Original Code written by: Thomakos, D. based on Shumway and Stoffer.
##------------------------------------------------------------------------
## Edited by: Fotis Papailias (f.papailias@qmul.ac.uk)
##------------------------------------------------------------------------
## Notes: This method works fine for more complec ARFIMA (p, d, q) models.
##------------------------------------------------------------------------
## Reference:
##------------------------------------------------------------------------
## You can call the function, in general, as "whittle2(x)" where "x" is the
## underlying time series.
##------------------------------------------------------------------------

logLik.Whittle <- function(d, xeval, demean = TRUE)
{
    Q <- spec.pgram(xeval, demean = demean, plot = FALSE,
                    na.action = na.omit)
    P <- Q$spec
    f <- Q$freq
    M <- length(f)
    g <- 4*(sin(pi*f)^2)
    s2 <- mean((g^d)*P)
    M*log(s2) + d*sum(log(g)) + M
}

whittle2 <- function(x, demean = TRUE, tol = 1e-4)
{
    ## Estimate the long memory parameter of the demeaned series
    if(demean) x <- x - mean(x)
    f <- optimize(logLik.Whittle, interval = c(-0.5,0.5),
                  xeval = x, demean = FALSE, tol = tol)
    list(d = f$minimum, loglik = f$objective)
}

