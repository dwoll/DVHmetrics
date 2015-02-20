## LOESS loss function - AICc or GCV
loessLoss <- function(span, x, method=c("AICc", "GCV"), data) {
    stopifnot(inherits(x, "loess"))
    method <- match.arg(method)

    fit     <- update(x, span=span, data=data)
    n       <- fit$n
    tr      <- fit$trace.hat
    sigmaSq <- sum(residuals(fit)^2) / (fit$n-1)

    if(method == "AICc") {
        ## AICc
        ## Hurvich, CM, Simonoff, JS, and Tsai, CL. 1998.
        ## Smoothing parameter selection in nonparametric regression using
        ## an improved Akaike Information Criterion.
        ## Journal of the Royal Statistical Society B, 60, 271-293.
        log(sigmaSq) + 1 + 2 * (2*(tr+1)) / (n-tr-2)
    } else if(method == "GCV") {
        ## generalized cross-validation
        ## Takezawa, K. 2005.
        ## Introduction to Nonparametric Regression. Wiley.
        n*sigmaSq / (n-tr)^2
    }
}

## LOESS fit with minimal AICc or GCV
loessSelSpan <- function(x, data, spanLoHi=c(0.1, 0.9), method=c("AICc", "GCV")) {
    stopifnot(inherits(x, "loess"), length(spanLoHi) == 2L)
    method <- match.arg(method)

    ## span for optimal loess fit
    span <- optimize(loessLoss, interval=spanLoHi, x=x, method=method, data=data)$minimum
    update(x, span=span, data=data)
}
