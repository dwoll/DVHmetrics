## we need ... because getMetric() will also pass parameters
## intended for other functions through ...
getEUD <-
function(x, EUDa, EUDfd=NULL, EUDab=NULL, ...) {
    UseMethod("getEUD")
}

getEUD.DVHs <-
function(x, EUDa, EUDfd=NULL, EUDab=NULL, ...) {
    if(length(EUDa) > 1L) {
		warning(paste0("Will only use EUDa=", EUDa[1]))
		EUDa <- EUDa[1]
	}

    if(isTRUE(all.equal(EUDa, 0))) {
		warning("'EUDa' must not be zero")
		return(NA_real_)
	}

    ## special cases
	gEUD <- if(isTRUE(all.equal(EUDa, 1))) {
		getMetric(x, "DMEAN")$DMEAN
	} else if(is.infinite(EUDa) && (EUDa > 0)) {  # +Inf
		getMetric(x, "DMAX")$DMAX
	} else if(is.infinite(EUDa) && (EUDa < 0)) {  # -Inf
		getMetric(x, "DMIN")$DMIN
	} else {
        ## get differential DVH
        xD <- convertDVH(x, toType="differential", toDoseUnit="asis", perDose=FALSE)
    
        ## convert dose to EQD2 if possible
        volume <- xD$dvhDiff[ , "volume"]
        dose   <- if(!is.null(EUDfd) && !is.null(EUDab)) {
            getEQD2(D=xD$dvhDiff[ , "dose"], fd=EUDfd, ab=EUDab)$EQD2
        } else {
            xD$dvhDiff[ , "dose"]
        }
    
        ## numerically unstable with large dose (in cGy) and EUDa
        ## -> take logs log(volume) - log(xD$structVol) + EUDa * log(dose)
        ## and try to carry logs up to geud calculation, only then go
        ## back to original scale
        volDose <- (volume / xD$structVol) * (dose^EUDa)
        wtMean  <- sum(volDose[volume > 0], na.rm=TRUE)
        geud    <- wtMean^(1/EUDa)
    	if(!is.finite(geud)) {
    		warning("Numerical problems encountered, NA returned")
    		geud <- NA_real_
    	}
        
        geud
	}

    data.frame(EUD=gEUD,
               patID=x$patID,
               structure=x$structure,
               stringsAsFactors=FALSE)
}

getEUD.DVHLst <-
function(x, EUDa, EUDfd=NULL, EUDab=NULL, ...) {
    EUDl  <- Map(getEUD, x,
                 EUDa=list(EUDa), EUDfd=list(EUDfd), EUDab=list(EUDab))
    EUDdf <- do.call("rbind", EUDl)
    rownames(EUDdf) <- NULL
    EUDdf
}

getEUD.DVHLstLst <-
function(x, EUDa, EUDfd=NULL, EUDab=NULL, ...) {
    EUDl  <- Map(getEUD, x,
                 EUDa=list(EUDa), EUDfd=list(EUDfd), EUDab=list(EUDab))
    EUDdf <- do.call("rbind", EUDl)
    rownames(EUDdf) <- NULL
    EUDdf
}
