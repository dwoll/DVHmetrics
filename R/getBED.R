#####---------------------------------------------------------------------------
## implement recycling rule for function arguments
#####---------------------------------------------------------------------------

recycle <-
function(...) {
    dots <- list(...)
    maxL <- max(sapply(dots, length))
    lapply(dots, rep, length=maxL)
}

#####---------------------------------------------------------------------------
## BED, EQD2, isoeffective dose
#####---------------------------------------------------------------------------

getBED <-
function(D=NULL, fd=NULL, fn=NULL, ab=NULL) {
    stopifnot(!is.null(fd), !is.null(ab))

    argL <- suppressWarnings(recycle(D, fd, fn, ab))
    D    <- argL[[1]]
    fd   <- argL[[2]]
    fn   <- argL[[3]]
    ab   <- argL[[4]]

    keepAB <- ab > 0
    if(any(!keepAB)) {
		warning("'ab' must be > 0")
    }
        
    keepFD <- fd > 0
    if(any(!keepFD)) {
        warning("'fd' must be > 0")
    }

    if(is.null(D)) {
        if(is.null(fn)) {
            stop("Either 'D' or 'fn' must be specified")
        } else {
            fn <- as.integer(fn)
            keepFN <- fn > 0
            if(any(!keepFN)) {
                warning("'fn' must be an integer > 0")
            }

            keep <- keepAB & keepFD & keepFN
            D    <- fn*fd
        }
    } else {
        keepD <- D > 0
        if(any(!keepD)) {
        	warning("'D' must be > 0")
        }

        keep <- keepD & keepAB & keepFD
    }

    BED <- rep(NA_real_, times=length(D))
    BED[keep] <- D[keep] * (1 + (fd[keep]/ab[keep]))
    return(BED)
}

getEQD2 <-
function(D=NULL, fd=NULL, fn=NULL, ab=NULL) {
    stopifnot(!is.null(fd), !is.null(ab))

    argL <- suppressWarnings(recycle(D, fd, fn, ab))
    D    <- argL[[1]]
    fd   <- argL[[2]]
    fn   <- argL[[3]]
    ab   <- argL[[4]]

    keepAB <- ab > 0
    if(any(!keepAB)) {
    	warning("'ab' must be > 0")
    }
        
    keepFD <- fd > 0
    if(any(!keepFD)) {
	    warning("'fd' must be > 0")
    }

    if(is.null(D)) {
        if(is.null(fn)) {
            stop("Either 'D' or 'fn' must be specified")
        } else {
            fn <- as.integer(fn)
            keepFN <- fn > 0
            if(any(!keepFN)) {
                warning("'fn' must be an integer > 0")
            }

            keep <- keepAB & keepFD & keepFN
            D    <- fn*fd
        }
    } else {
        keepD <- D > 0
        if(any(!keepD)) {
        	warning("'D' must be > 0")
        }

        keep <- keepD & keepAB & keepFD
    }

    BED <- D * (1 + (fd/ab))

    ## EQD2 <- D * (fd + ab) / (2 + ab)
    EQD2 <- rep(NA_real_, times=length(D))
    EQD2[keep] <- BED[keep] / (1 + (2/ab[keep]))

    return(EQD2)
}

getIsoEffD <-
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
