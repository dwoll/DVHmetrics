#####---------------------------------------------------------------------------
## EQD2 linear quadratic model
#####---------------------------------------------------------------------------

getEQD2 <-
function(D=NULL, fd=NULL, fn=NULL, ab=NULL) {
    UseMethod("getEQD2")
}

getEQD2.default <-
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

    data.frame(EQD2=EQD2,
               fractDose=fd,
               ab=ab)
}

getEQD2.DVHs <-
function(D=NULL, fd=NULL, fn=NULL, ab=NULL) {
    stopifnot(!is.null(D), !is.null(fn))

    Dmax <- getMetric(D, "DMAX")$DMAX
    fd   <- Dmax/fn
    D$dvh[ , "dose"] <- getEQD2(D=D$dvh[ , "dose"], fd=fd, ab=ab)$EQD2
    D
}

getEQD2.DVHLst <-
function(D=NULL, fd=NULL, fn=NULL, ab=NULL) {
    stopifnot(!is.null(D))

    EQD2l <- Map(getEQD2, D, fd=list(fd), fn=list(fn), ab=list(ab))
    class(EQD2l) <- "DVHLst"
    attr(EQD2l, which="byPat") <- TRUE
    EQD2l
}

getEQD2.DVHLstLst <-
function(D=NULL, fd=NULL, fn=NULL, ab=NULL) {
    stopifnot(!is.null(D))

    EQD2ll <- Map(getEQD2, D, fd=list(fd), fn=list(fn), ab=list(ab))
    class(EQD2ll) <- "DVHLstLst"
    attr(EQD2ll, which="byPat") <- TRUE
    EQD2ll
}
