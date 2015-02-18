## parse metric character strings into components and ensure validity
## optionally convert dose unit and volume unit
parseMetric <- function(x, doseUnit=NULL, volUnit=NULL) {
    ## remove whitespace and convert to upper case
    x <- toupper(removeWS(x))
    
    ## regular expression for components of a metric character string
    ## allow for V10%_CC or D10cc_% pattern for returning absolute volumes / relative doses
    pattern <- "^([DV])([.[:digit:]]+|MAX|MIN|MEAN|RX|SD)([%]|GY|CGY|CC)*(_[%]|_GY|_CGY|_CC)*.*"
    
    ## extract components
    DV      <- sub(pattern, "\\1", x)   # get a volume or a dose?
    valRef  <- sub(pattern, "\\2", x)   # given value at which to evaluate volume or dose
    unitRef <- sub(pattern, "\\3", x)   # measurement unit for given volume or dose
    unitDV_ <- sub(pattern, "\\4", x)   # measurement unit for output volume or dose

    ## remove _ from unitDV_
    unitDV <- ifelse(unitDV_ != "", sub("_", "", unitDV_), NA_character_)

    ## special dose cases: DMEAN, DSD, DMIN, DMAX, DRX
    specDose  <- valRef %in% c("MAX", "MIN", "MEAN", "RX", "SD")
    valRefNum <- ifelse((DV == "D") & specDose, NA_real_,      suppressWarnings(as.numeric(valRef)))
    unitRef   <- ifelse((DV == "D") & specDose, NA_character_, unitRef)

    ## convert absolute dose units if requested
    if(!is.null(doseUnit)) {
        doseUnit <- toupper(removeWS(doseUnit))

        ## output is dose -> reference is volume, just set unitDV
        ## don't convert relative volume
        idxD   <- (DV == "D") & (unitDV != "%")
        unitDV <- ifelse(idxD, doseUnit, unitDV)

        ## output is volume -> reference is dose
        ## don't convert relative dose
        idxV      <- (DV == "V") & (unitRef != "%")
        valRefNum <- ifelse(idxV,
                            valRefNum*suppressWarnings(getConvFac(paste0(unitRef, "2", doseUnit))),
                            valRefNum)

        valRef  <- ifelse(idxV, as.character(valRefNum), valRef)
        unitRef <- ifelse(idxV, doseUnit, unitRef)
    }
    
    ## convert absolute volume units if requested
    if(!is.null(volUnit)) {
        volUnit <- toupper(removeWS(volUnit))

        ## output is dose -> reference is volume
        ## don't convert relative volume
        ## consider special dose cases
        idxD      <- (DV == "D") & (unitRef != "%") & !specDose
        valRefNum <- ifelse(idxD,
                            valRefNum*suppressWarnings(getConvFac(paste0(unitRef, "2", volUnit))),
                            valRefNum)

        valRef  <- ifelse(idxD, as.character(valRefNum), valRef)
        unitRef <- ifelse(idxD, volUnit, unitRef)

        ## output is volume -> reference is dose, just set unitDV
        ## don't convert relative dose
        idxV   <- (DV == "V") & (unitDV != "%")
        unitDV <- ifelse(idxV, volUnit, unitDV)
    }

    ## canonical metric string
    unitDVStr  <- ifelse(is.na(unitDV),  "", paste0("_", unitDV))
    unitRefStr <- ifelse(is.na(unitRef), "", unitRef)
    metric     <- paste0(DV, valRef, unitRefStr, unitDVStr)

    ## check validity
    ## consider special cases DMEAN, DSD, DMIN, DMAX, DRX
    ## V -> %, Gy, cGy
    ## D -> %, CC
    validPattern <- grepl(pattern, x)
    validUnitRef <- ((DV == "D") &
                     ((unitRef %in% c("%", "CC")) |
                      (valRef  %in% c("MAX", "MIN", "MEAN", "RX", "SD")))) |
                    ((DV == "V") & (unitRef %in% c("%", "GY", "CGY")))
    validUnitDV  <- is.na(unitDV) |
                    ((DV == "D") & (unitDV  %in% c("%", "GY", "CGY"))) |
                    ((DV == "V") & (unitDV  %in% c("%", "CC")))

    valid <- validPattern & validUnitRef & validUnitDV
    if(!all(valid)) {
        warning(c("Pattern ", paste(x[!valid], collapse=", "), " is invalid"))
    }
 
    return(data.frame(metric=metric, valid=valid,
                      DV=DV, unitDV=unitDV,
                      valRef=valRef, valRefNum=valRefNum, unitRef=unitRef,
                      stringsAsFactors=FALSE))
}

#####---------------------------------------------------------------------------
## S3 generic function
getMetric <-
function(x, metric, patID, structure,
     sortBy=c("none", "observed", "patID", "structure", "metric"),
    splitBy=c("none", "patID", "structure", "metric"),
     interp=c("linear", "spline"), ...) {
    UseMethod("getMetric")
}

## get metrics for one DVH
## patID, structure, sortBy, splitBy are ignored
getMetric.DVHs <-
function(x, metric, patID, structure,
     sortBy=c("none", "observed", "patID", "structure", "metric"),
    splitBy=c("none", "patID", "structure", "metric"),
     interp=c("linear", "spline"), ...) {

    interp <- match.arg(interp)

    ## D*%  = min radiation of most extremely radiated *% volume
    ## D*cc = min radiation of most extremely radiated *cc volume
    getDose <- function(val, type="relative", unitRef, unitDV) {
        vol <- if(type == "relative") {
            x$dvh[ , "volumeRel"]
        } else if(type == "absolute") {
            x$dvh[ , "volume"]
        }

        ## warn if beyond smallest/largest available volume
        ## like in D0.001cc, D100000cc
        if(isTRUE((val < min(vol)) || (val > max(vol)))) {
            warning(c("Requested ", type, " volume\n", val,
                      " is outside the DVH interval [",
                      formatC(min(vol), format="f", digits=3), ",",
                      formatC(max(vol), format="f", digits=3), "]\n",
                      "NA returned\n"))
        }

        dose <- if(!is.na(unitDV) && (unitDV == "%")) {   # return relative dose
            x$dvh[ , "doseRel"]
        } else {                         # return absolute dose (default)
            cf <- if(!is.na(unitDV)) {
                getConvFac(paste0(x$doseUnit, "2", unitDV))
            } else {
                1
            }
            cf * x$dvh[ , "dose"]
        }

        if(interp == "linear") {  # rule=1 -> no interpolation beyond bounds
            res <- try(approx(vol, dose, val, method="linear", rule=1)$y)
            if(!inherits(res, "try-error")) {
                res
            } else {
                NA_real_
            }
        } else if(interp == "spline") {  # does interpolation beyond bounds
            spline(vol, dose, val, method="hyman")$y
        }
    }

    ## V*%  = max volume with radiation >= *% of prescribed dose
    ## V*Gy = max volume with radiation >= * Gy / cGy
    getVolume <- function(val, type="relative", unitRef, unitDV) {
        if(type == "relative") {         # given dose is relative
            dose <- x$dvh[ , "doseRel"]
            val  <- val/100

            ## check if max dose is smaller than requested % of prescribed dose
            if(val > 1) {
                warning("max dose is less than requested dose")
                return(0)
            }
        } else if(type == "absolute") {  # given dose is absolute
            cf   <- getConvFac(paste0(x$doseUnit, "2", unitRef))
            dose <- cf * x$dvh[ , "dose"]

            ## check if max dose is smaller than requested % of prescribed dose
            if(val > (cf * x$doseMax)) {
                warning("max dose is less than requested dose")
                return(0)
            }
        }

        vol <- if(!is.na(unitDV) && (unitDV == "CC")) {    # return absolute volume
            x$dvh[ , "volume"]
        } else {                         # return relative volume (default)
            x$dvh[ , "volumeRel"]
        }

        if(interp == "linear") {  # rule=1 -> no interpolation beyond bounds
            res <- try(approx(dose, vol, val, method="linear", rule=1)$y)
            if(!inherits(res, "try-error")) {
                res
            } else {
                NA_real_
            }
        } else if(interp == "spline") {  # does interpolation beyond bounds
            spline(dose, vol, val, method="hyman")$y
        }
    }

    ## get value for 1 parsed metric
    getVal <- function(pm) {
        if(!pm$valid) {
            return(NA_real_)
        } else if(pm$DV == "D") {              # report a dose
            if(pm$valRef %in% c("MAX", "MIN", "MEAN", "RX", "SD")) {
                ## mean from differential DVH from
                ## integration of dose * (spline fit derivative)
                ## this breaks for very fine-grained DVHs
                # xx  <- seq(0, max(dose), length.out=1000)
                # spl <- smooth.spline(dose, 1 - x$dvh[ , "volumeRel"])
                # dy  <- predict(spl, xx, deriv=1)$y
                # meanSpl <- tryCatch(integrate(function(y) {
                #    y*predict(spl, y, deriv=1)$y/100 }, 0, max(dose))$value,
                #    error=function(e) return(NA))
                cf <- if(!is.na(pm$unitDV)) {
                    getConvFac(paste0(x$doseUnit, "2", pm$unitDV))
                } else {
                    1
                }

                cf * switch(pm$valRef,
            		MIN = x$doseMin,
        			MAX = x$doseMax,
    				MEAN= x$doseAvg,       # meanMid, meanSpl
    				RX  = x$doseRx,
                    SD  = x$doseSD)
            } else if(pm$unitRef == "%") {
                getDose(pm$valRefNum,   type="relative",
                        unitRef=pm$unitRef, unitDV=pm$unitDV)
            } else if(pm$unitRef == "CC") {
                getDose(pm$valRefNum,   type="absolute",
                        unitRef=pm$unitRef, unitDV=pm$unitDV)
            }
        } else if(pm$DV == "V") {          # report a volume
            if(pm$unitRef == "%") {
                getVolume(pm$valRefNum, type="relative",
                          unitRef=pm$unitRef, unitDV=pm$unitDV)
            } else if(pm$unitRef %in% c("GY", "CGY")) {
                getVolume(pm$valRefNum, type="absolute",
                          unitRef=pm$unitRef, unitDV=pm$unitDV)
            }
        }
    }
    
    ## parse metric strings into lists of components
    pm  <- parseMetric(metric)
    pmL <- split(pm, f=seq_len(nrow(pm)))
    pmL <- setNames(pmL, vapply(pmL, function(y) y$metric, character(1)))

    ## calculate results for metrics
    Map(getVal, pmL)
}

## getMetric.DVHLst(): return list of separate data frames
## one for each patient, structure, or metric
## x is a list of DVH objects - 1 for each structure or 1 for each patient
getMetric.DVHLst <-
function(x, metric, patID, structure,
     sortBy=c("none", "observed", "patID", "structure", "metric"),
    splitBy=c("none", "patID", "structure", "metric"),
     interp=c("linear", "spline"), ...) {

    if(!missing(patID)) {
        ## filter by patID
        p <- trimWS(patID, side="both")
        x <- Filter(function(y) any(grepl(paste(p, collapse="|"), y$patID, ...)), x)
        if(length(x) < 1L) { stop("selected patient not found") }
    }

    if(!missing(structure)) {
        ## filter by structure
        s <- trimWS(structure, side="both")
        x <- Filter(function(y) any(grepl(paste(s, collapse="|"), y$structure, ...)), x)
        if(length(x) < 1L) { stop("selected structure not found") }
    }

    res    <- Map(getMetric, x, metric=list(metric), interp=list(interp))
    resDFL <- melt(res, value.name="observed")
    names(resDFL)[names(resDFL) == "L1"] <- "structure"
    names(resDFL)[names(resDFL) == "L2"] <- "metric"

    ## sort if requested
    if(!("none" %in% sortBy)) {
        idx    <- do.call("order", lapply(sortBy, function(y) with(resDFL, get(y))))
        resDFL <- resDFL[idx, ]
    }

    ## split if requested
    if(!("none" %in% splitBy)) {
        splRes <- split(resDFL, lapply(splitBy, function(y) with(resDFL, get(y))))
        ## only keep non-empty components
        resDFL <- Filter(function(y) nrow(y) > 0, splRes)
    }

    ## if result is list of length 1 -> convert to data frame
    if((!is.data.frame(resDFL)) && (length(resDFL) == 1)) {
        resDFL <- resDFL[[1]]
    }

    return(resDFL)
}

## getMetric.DVHLstLst(): return list of separate data frames
## one for each structure, or one for each metric
## x is from many DVH-files -> list of DVH lists
getMetric.DVHLstLst <-
function(x, metric, patID, structure,
     sortBy=c("none", "observed", "patID", "structure", "metric"),
    splitBy=c("none", "patID", "structure", "metric"),
     interp=c("linear", "spline"), ...) {

    ## re-organize x into by-patient form if necessary
    isByPat <- attributes(x)$byPat
    if(is.null(isByPat) || !isByPat) {
        x <- reorgByPat(x, byPat=TRUE)
    }

    ## if patIDs are selected, filter those
    if(!missing(patID)) {
        p <- trimWS(patID, side="both")
        x <- Filter(function(y) any(grepl(paste(p, collapse="|"), y[[1]]$patID, ...)), x)
        if(length(x) < 1L) { stop("No selected patient found") }
    }

    struct <- if(missing(structure)) {
        NULL
    } else {
        structure
    }

    ## DVH list y is from 1 DVH-file: list of DVHs - many structures for one ID
    collectMetrics <- function(y, metric, structure=NULL) {
        ## if structures are selected, filter those
        if(!is.null(structure)) {
            s <- trimWS(structure, side="both")
            y <- Filter(function(z) any(grepl(paste(s, collapse="|"), z$structure, ...)), y)
            if(length(y) < 1L) { stop("No selected structure found") }
        }

        Map(getMetric, y, metric=list(metric), interp=list(interp))
    }

    res    <- Map(collectMetrics, x, metric=list(metric), structure=list(struct))
    resDFL <- melt(res, value.name="observed")
    names(resDFL)[names(resDFL) == "L1"] <- "patID"
    names(resDFL)[names(resDFL) == "L2"] <- "structure"
    names(resDFL)[names(resDFL) == "L3"] <- "metric"

    ## sort if requested
    if(!("none" %in% sortBy)) {
        idx    <- do.call("order", lapply(sortBy, function(y) with(resDFL, get(y))))
        resDFL <- resDFL[idx, ]
    }

    ## split if requested
    if(!("none" %in% splitBy)) {
        splRes <- split(resDFL, lapply(splitBy, function(y) with(resDFL, get(y))))
        ## only keep non-empty components
        resDFL <- Filter(function(y) nrow(y) > 0, splRes)
    }

    ## if result is list of length 1 -> convert to data frame
    if((!is.data.frame(resDFL)) && (length(resDFL) == 1)) {
        resDFL <- resDFL[[1]]
    }

    return(resDFL)
}
