## we need ... because getMetric() will also pass parameters
## intended for other functions through ...
getEUD <-
function(x, EUDa, EUDfn=NULL, EUDab=NULL, ...) {
    UseMethod("getEUD")
}

getEUD.DVHs <-
function(x, EUDa, EUDfn=NULL, EUDab=NULL, ...) {
    if(length(EUDa) > 1) {
		warning(paste0("Will only use EUDa=", EUDa[1]))
		EUDa <- EUDa[1]
	}

    if(isTRUE(all.equal(EUDa, 0))) {
		warning("'EUDa' must not be zero")
		return(NA_real_)
	}

    ## special cases
	if(isTRUE(all.equal(EUDa, 1))) {
		return(getMetric(x, "DMEAN")$DMEAN)
	} else if(is.infinite(EUDa) && (EUDa > 0)) {  # +Inf
		return(getMetric(x, "DMAX")$DMAX)
	} else if(is.infinite(EUDa) && (EUDa < 0)) {  # -Inf
		return(getMetric(x, "DMIN")$DMIN)
	}

    ## get differential DVH
    xD <- convertDVH(x, toType="differential", toDoseUnit="asis", perDose=FALSE)
    
    ## convert dose to EQD2 if possible
    volume <- xD$dvh[ , "volume"]
    dose   <- if(!is.null(EUDfn) && !is.null(EUDab)) {
        Dmax <- getMetric(x, "DMAX")$DMAX
        fd   <- Dmax/EUDfn
        getEQD2(xD$dvh[ , "dose"], D=Dmax, fd=fd, ab=EUDab)
    } else {
        xD$dvh[ , "dose"]
    }

    volDose <- volume*dose^EUDa / xD$structVol
    wtMean  <- sum(volDose[volume > 0], na.rm=TRUE)
    gEUD    <- wtMean^(1/EUDa)
	if(!is.finite(gEUD)) {
		warning("Numerical problems encountered, NA returned")
		gEUD <- NA_real_
	}

    return(gEUD) 
}

getEUD.DVHLst <-
function(x, EUDa, EUDfn=NULL, EUDab=NULL, ...) {
    unlist(Map(getEUD, x,
               EUDa=list(EUDa), EUDfn=list(EUDfn), EUDab=list(EUDab)), recursive=FALSE)
}

getEUD.DVHLstLst <-
function(x, EUDa, EUDfn=NULL, EUDab=NULL, ...) {
    unlist(Map(getEUD, x,
               EUDa=list(EUDa), EUDfn=list(EUDfn), EUDab=list(EUDab)), recursive=FALSE)
}