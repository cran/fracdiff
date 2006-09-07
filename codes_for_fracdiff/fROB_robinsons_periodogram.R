## Robinson's Log Periodogram Estimation Method
##------------------------------------------------------------------------
## Original Code written by: Fotis Papailias
## based on and using parts of the codes written by Reisen, V. for the GPH
## and Sperio estimators of the 'fracdiff' package.
##------------------------------------------------------------------------
## Notes: This method introduced by Robinson (1995) and it is
## theoritically, more improved than GPH. However the proper choice of the
## bandwidth exponent and the other parameters that requires can cause problems.
##------------------------------------------------------------------------
## References:
##
## [1]
## Robinson, P. M. (1995): “Log-Peridogram Regression of Time Series with
## Long Range Dependence”, Annals of Statistics, 23 (3), 1048-1072.
##
## [2]
## Reisen, V. A. and S. Lopes (1999): “Some Simulations and Applications of Forecasting Long-
## Memory Time Sereis Models”, Journal of Statistical Planning and Inference, 80, 269-287.
##
## [3]
## Reisen, V. A., Abraham B. and S. Lopes (2001): “Estimation of Parameters in ARFIMA
## Processes: A Simulation Study”, Communication in Statistics: Simulation and Computation, 30
## (4), 787-803.
##
## [4]
## Reisen, V. A., Abraham, B. and E. M. Toscano (2000): “Parametric and Semi Parametric
## Estimations of Stationary Univariate ARFIMA Models”,  Brazilian Journal of Probability and
## Statistics, 14, 185-206.
##------------------------------------------------------------------------
## You can call the function, in general, as "fROB(x, d, el, alpha, tau)"
## where "x" is the underlying time series, "el" the starting lag to be used and
## "alpha" and "tau" define the equation for the last lag to be used.
##------------------------------------------------------------------------
## ATTENTION!!
##------------------------------------------------------------------------
## Because of the complexity of the "lag equation" study [2] or [3] or [4]
## from the above references and observe how the authors define it.
## Then, you would set the parameters very easy.
##------------------------------------------------------------------------


fROB <- function(x, d, el, alpha, tau)
{
    x <- as.numeric(na.fail(as.ts(x1)))
    if (any(is.na(x))) stop("NAs in x")
    if (NCOL(x) > 1) stop("only implemented for univariate time series")
    n <- length(x)
    g <-
        if (d >= 0 & d <= 0.25) {
            alpha*(n^((2*tau)/((2*tau)+1)))
        } else if (d > 0.25 & d < 0.5) {
            alpha*(n^(tau/(tau+1-(2*d))))
        } else 0
    j <- el:g
    kk <- 1:(n-1)
    w <- 2*pi*j/n
    mx <- mean(x)
    x <- x - mx
    ##--- begin{FIXME} --- use fft() [or spec.pgram()] instead of the following:
    var.x <- sum(x^2)/n ## not /(n-1)
    cov.x <- numeric(n-1)
    for (k in kk)
        cov.x[k] <- sum(x[1:(n-k)] * x[(1+k):n]) / n
    periodogram <- numeric(g-el)
    z <- g-(el-1)
    for (i in 1:z)
        periodogram[i] <- var.x + 2*sum(cov.x * cos(w[i]*kk))
    y.reg <- log(periodogram / (2*pi))
    x.reg <- 2*log(2*sin(w/2))
    ##--- end{FIXME} ---
    fit <- lm(y.reg ~ x.reg) ## MM: "FIXME" - use lm.fit()
    dROB <- coef(fit)[2]
    names(dROB) <- NULL
    x.r2 <- sum((x.reg - mean(x.reg))^2)
    var.d <- pi^2 / (6*x.r2)
    var.reg <- sum(resid(fit)^2) / ((g - 1) * x.r2)
    list(d = -dROB, sd.as = sqrt(var.d), sd.reg = sqrt(var.reg))
}
