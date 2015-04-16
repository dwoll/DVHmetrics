#####---------------------------------------------------------------------------
## isoeffective dose linear quadratic model
#####---------------------------------------------------------------------------

getIsoEffD <-
function(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, fn1=NULL, fn2=NULL, ab=NULL) {
    UseMethod("getIsoEffD")
}

getIsoEffD.default <-
function(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, fn1=NULL, fn2=NULL, ab=NULL) {
    stopifnot(!is.null(fd1), !is.null(ab))

    argL <- suppressWarnings(recycle(D1, D2, fd1, fd2, fn1, fn2, ab))
    maxL <- max(sapply(argL, length))
    D1   <- argL[[1]]
    D2   <- argL[[2]]
    fd1  <- argL[[3]]
    fd2  <- argL[[4]]
    fn1  <- argL[[5]]
    fn2  <- argL[[6]]
    ab   <- argL[[7]]

    keepFD1 <- fd1 > 0
    if(any(!keepFD1)) {
        warning("'fd1' must be > 0")
    }

    keepAB <- ab > 0
    if(any(!keepAB)) {
		warning("'ab' must be > 0")
    }

    if(is.null(D1)) {
        if(is.null(fn1)) {
            stop("Without D1, need fd1 and fn1")
        }

        fn1 <- as.integer(fn1)
        keepFN1 <- fn1 > 0
        if(any(!keepFN1)) {
            warning("'fn1' must be an integer > 0")
        }

        D1 <- fd1*fn1
    }
    
    keepD1 <- D1 > 0
    if(any(!keepD1)) {
    	warning("'D1' must be > 0")
    }

    ## cases
    if(is.null(D2) && is.null(fn2) && !is.null(fd2)) {
        ## D1 (or fd1, fn1) are given, D2 and fn2 are missing
        ## -> need fd1, fd2, ab -> calculate D2 (special case: EQD2 with fd2=2)
        keep <- keepD1 & keepFD1 & keepAB
        D2   <- rep(NA_real_, times=maxL)
        D2[keep] <- D1[keep] * (fd1[keep] + ab[keep]) / (fd2[keep] + ab[keep])
        D2
    } else if(is.null(fd2)) {
        ## D1 (or fd1, fn1) and D2 are given, fd2 is missing
        ## -> need fd1, ab -> calculate fd2
        if(is.null(D2)) {
            if(is.null(fn2)) {
                stop("Without D2, need fd2 and fn2")
            }
    
            fn2 <- as.integer(fn2)
            keepFN2 <- fn2 > 0
            if(any(!keepFN2)) {
                warning("'fn2' must be an integer > 0")
            }
    
            D2 <- fd2*fn2
        }
        
        keepD2 <- D2 > 0
        if(any(!keepD2)) {
        	warning("'D2' must be > 0")
        }

        keep <- keepD1 & keepD2 & keepFD1 & keepAB
        fd2  <- rep(NA_real_, times=maxL)
        fd2[keep] <- (D1[keep] / D2[keep]) * (fd1[keep] + ab[keep]) - ab[keep]
        fd2
    } else {
        rep(NA_real_, times=maxL)
    }
}

getIsoEffD.DVHs <-
function(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, fn1=NULL, fn2=NULL, ab=NULL) {
    stopifnot(!is.null(D1), is.null(D2), is.null(fd1), !is.null(fd2), !is.null(fn1), is.null(fn2))

    Dmax <- getMetric(D1, "DMAX")$DMAX
    fd1  <- Dmax/fn1
    D1$dvh[ , "dose"] <- getIsoEffD(D1=D1$dvh[ , "dose"], D2=D2, fd1=fd1, fd2=fd2,
                                    fn1=fn1, fn2=fn2, ab=ab)
    D1
}

getIsoEffD.DVHLst <-
function(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, fn1=NULL, fn2=NULL, ab=NULL) {
    stopifnot(!is.null(D1))

    IsoEDl <- Map(getIsoEffD, D1, fd1=list(fd1), fd2=list(fd2),
                  fn1=list(fn1), fn2=list(fn2), ab=list(ab))
    class(IsoEDl) <- "DVHLst"
    attr(IsoEDl, which="byPat") <- TRUE
    IsoEDl
}

getIsoEffD.DVHLstLst <-
function(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, fn1=NULL, fn2=NULL, ab=NULL) {
    stopifnot(!is.null(D1))

    IsoEDll <- Map(getIsoEffD, D1, fd1=list(fd1), fd2=list(fd2),
                   fn1=list(fn1), fn2=list(fn2), ab=list(ab))
    class(IsoEDll) <- "DVHLstLst"
    attr(IsoEDll, which="byPat") <- TRUE
    IsoEDll
}
