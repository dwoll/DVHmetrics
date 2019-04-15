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

        volume <- volume[volume >= 0]
        dose   <- dose[volume >= 0]

        ## numerically unstable with large dose (in cGy) and EUDa
        # volDose <- (volume / xD$structVol) * (dose^EUDa)
        # wtMean  <- sum(volDose[volume > 0], na.rm=TRUE)
        # wtMean  <- weighted.mean(dose^EUDa, w=volume, na.rm=TRUE)
        # geud    <- wtMean^(1/EUDa)
        ## rescale dose first, then scale back weighted mean
        ## https://stackoverflow.com/a/47296640/484139
        maxVal <- if(EUDa > 0) {
            max((volume/xD$structVol)^(1/EUDa) * dose)
        } else if(EUDa < 0) { # avoid 1 / 0^EUDa problem
            max(dose)
        } else {
        	stop("EUDa must not be 0")
        }

        doseScale <- dose / maxVal
        wtMean    <- weighted.mean(doseScale^EUDa, w=volume, na.rm=TRUE)
        geud      <- maxVal * wtMean^(1/EUDa)
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
