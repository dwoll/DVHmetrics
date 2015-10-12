## convert DVH to cumulative/differential DVH
convertDVH <-
function(x, toType=c("asis", "cumulative", "differential"),
         toDoseUnit=c("asis", "GY", "CGY"),
         interp=c("asis", "linear", "spline", "ksmooth", "smoothSpl"),
         nodes=NULL, perDose=TRUE) {
    UseMethod("convertDVH")
}

## for DVH matrix itself
convertDVH.matrix <-
function(x, toType=c("asis", "cumulative", "differential"),
         toDoseUnit=c("asis", "GY", "CGY"),
         interp=c("asis", "linear", "spline", "ksmooth", "smoothSpl"),
         nodes=NULL, perDose=TRUE) {
    toType     <- match.arg(toType)
    toDoseUnit <- match.arg(toDoseUnit)
    interp     <- match.arg(interp)
    
    ## matrix may include duplicate rows
    x <- unique(x)

    ## when dose entries are duplicated with different volumes (cumulative DVH),
    ## pick row with max volume
    ## split matrix into rows with unique dose
    ## preserve matrix class and keep col names with split.data.frame()
    xL <- if(!all(is.na(x[ , "dose"]))) {
        split.data.frame(x, x[ , "dose"])
    } else {
        split.data.frame(x, x[ , "doseRel"])
    }
    
    ## function to keep row with max volume
    keepMaxVol <- function(xSub) {
        idx <- if(!all(is.na(xSub[ , "volume"]))) {
            xSub[ , "volume"] == max(xSub[ , "volume"])
        } else{
            xSub[ , "volumeRel"] == max(xSub[ , "volumeRel"])
        }
        
        xSub[idx, ]
    }

    ## apply keepMaxVol and re-bind to matrix
    x <- do.call("rbind", lapply(xL, keepMaxVol))
    rownames(x) <- NULL

    ## number of DVH points -> determine nodes for interpolation
    N     <- nrow(x)
    nodes <- max(nodes, N)
    
    ## absolute and relative doses/volumes
    dose      <- x[ , "dose"]
    doseRel   <- x[ , "doseRel"]
    volume    <- x[ , "volume"]
    volumeRel <- x[ , "volumeRel"]
    
    ## convert dose unit
    doseConv <- if(toDoseUnit == "asis") {
        dose              ## nothing to do
    } else if(toDoseUnit == "GY") {
        dose / 100        ## cGy to Gy
    } else if(toDoseUnit == "CGY") {
        dose * 100        ## Gy to cGy
    }

    ## convert DVH type
    if(toType == "asis") {
        ## nothing to convert
        ## interpolate?
        if(!(interp %in% c("asis", "linear"))) {
            ## if differential && perDose == FALSE
            ## normalize -> interpolate -> re-normalize 
            ## check if volume is already sorted -> cumulative DVH
            volumeSel <- if(!any(is.na(volumeRel))) {
                volumeRel
            } else {
                volume
            }
    
            DVHtype <- if(isTRUE(all.equal(volumeSel, sort(volumeSel, decreasing=TRUE)))) {
                "cumulative"
            } else {
                "differential"
            }
            
            if((DVHtype == "differential") && !perDose) {
                ## if perDose == FALSE: normalize -> interpolate -> re-normalize
                binW         <- diff(c(-doseConv[1], doseConv))
                volumeNew    <- -diff(volume)    / binW
                volumeRelNew <- -diff(volumeRel) / binW
            }
            
            sm <- switch(interp,
                         ksmooth=getKSmooth(    doseConv, doseRel, volume, volumeRel, nodes=nodes),
                         smoothSpl=getSmoothSpl(doseConv, doseRel, volume, volumeRel, nodes=nodes),
                         getInterpSpl(          doseConv, doseRel, volume, volumeRel, nodes=nodes)) # default

            doseNew    <- sm[ , "dose"]
            doseRelNew <- sm[ , "doseRel"]

            if((DVHtype == "differential") && !perDose) {
                ## re-normalize to non-per-dose volume
                ## start bin at 0 because no out-of-range interpolation
                binW         <- diff(c(0, doseConv))
                volumeNew    <- sm[ , "volume"]    * binW
                volumeRelNew <- sm[ , "volumeRel"] * binW

            } else {
                volumeNew    <- sm[ , "volume"]
                volumeRelNew <- sm[ , "volumeRel"]
            }
        } else {
            doseNew      <- doseConv
            doseRelNew   <- doseRel
            volumeNew    <- volume
            volumeRelNew <- volumeRel
        }
    } else if(toType == "cumulative") {
        ## convert from differential to cumulative DVH
        ## interpolate?
        if(!(interp %in% c("asis", "linear"))) {
            if(!perDose) {
                ## if perDose == FALSE: normalize -> interpolate -> re-normalize
                binW         <- diff(c(-doseConv[1], doseConv))
                volumeNew    <- -diff(volume)    / binW
                volumeRelNew <- -diff(volumeRel) / binW
            }

            sm <- switch(interp,
                         ksmooth=getKSmooth(    doseConv, doseRel, volume, volumeRel, nodes=nodes),
                         smoothSpl=getSmoothSpl(doseConv, doseRel, volume, volumeRel, nodes=nodes),
                         getInterpSpl(          doseConv, doseRel, volume, volumeRel, nodes=nodes)) # default

            doseConv <- sm[ , "dose"]
            doseRel  <- sm[ , "doseRel"]
            
            if(!perDose) {
                ## re-normalize to non-per-dose volume
                ## start bin at 0 because no out-of-range interpolation
                binW      <- diff(c(0, doseConv))
                volume    <- sm[ , "volume"]    * binW
                volumeRel <- sm[ , "volumeRel"] * binW
            } else {
                volume    <- sm[ , "volume"]
                volumeRel <- sm[ , "volumeRel"]
            }
        }

        ## dose category half-widths
        doseCatHW    <- diff(doseConv)/2
        doseRelCatHW <- diff(doseRel)/2

        ## dose category mid-points, starting at 0
        doseMidPt    <- c(0, doseConv[-N] + doseCatHW)
        doseRelMidPt <- c(0, doseRel[-N]  + doseRelCatHW)
        
        ## add one more category beyond max dose
        doseNew    <- c(doseMidPt,    doseConv[N] + doseCatHW[N-1])
        doseRelNew <- c(doseRelMidPt, doseRel[N]  + doseRelCatHW[N-1])
        
        if(perDose) {
            ## differential DVH -> volume is per Gy -> mult with bin-width
            binW         <- diff(c(-doseConv[1], doseConv))
            volumeBin    <- volume*binW
            volumeRelBin <- volumeRel*binW
        } else {
            volumeBin    <- volume
            volumeRelBin <- volumeRel
        }

        ## volume "survival" curve starting at max volume
        volumeNew    <- sum(volumeBin)    - c(0, cumsum(volumeBin))
        volumeRelNew <- sum(volumeRelBin) - c(0, cumsum(volumeRelBin))
    } else if(toType == "differential") {
        ## convert from cumulative to differential DVH
        ## dose category half-widths
        doseCatHW    <- diff(doseConv)/2
        doseRelCatHW <- diff(doseRel)/2

        doseNew      <- doseConv[-length(doseConv)] + doseCatHW
        doseRelNew   <-  doseRel[-length(doseRel)]  + doseRelCatHW

        ## differential DVH -> volume is per Gy -> divide by bin-width
        binW <- 2*doseCatHW
        volumeNew <- if(perDose) {
            -diff(volume) / binW
        } else {
            -diff(volume)
        }

        volumeRelNew <- if(perDose) {
            -diff(volumeRel) / binW
        } else {
            -diff(volumeRel)
        }

        ## interpolate?
        if(!(interp %in% c("asis", "linear"))) {
            ## if perDose == FALSE: normalize -> interpolate -> re-normalize
            if(!perDose) {
                volumeNew    <- -diff(volume)    / binW
                volumeRelNew <- -diff(volumeRel) / binW
            }

            ## interpolate
            sm <- switch(interp,
                  ksmooth=getKSmooth(  doseNew, doseRelNew, volumeNew, volumeRelNew, nodes=nodes),
                smoothSpl=getSmoothSpl(doseNew, doseRelNew, volumeNew, volumeRelNew, nodes=nodes),
                getInterpSpl(          doseNew, doseRelNew, volumeNew, volumeRelNew, nodes=nodes)) # default

            doseNew    <- sm[ , "dose"]
            doseRelNew <- sm[ , "doseRel"]

            if(!perDose) {
                ## re-normalize to non-per-dose volume
                ## start bin at 0 because no out-of-range interpolation
                binW         <- diff(c(0, doseNew))
                volumeNew    <- sm[ , "volume"]   *binW
                volumeRelNew <- sm[ , "volumeRel"]*binW
            } else {
                volumeNew    <- sm[ , "volume"]
                volumeRelNew <- sm[ , "volumeRel"]
            }
        }
    }

    cbind(dose=doseNew, doseRel=doseRelNew,
          volume=volumeNew, volumeRel=volumeRelNew)
}

## for 1 DVH object
convertDVH.DVHs <-
function(x, toType=c("asis", "cumulative", "differential"),
         toDoseUnit=c("asis", "GY", "CGY"),
         interp=c("asis", "linear", "spline", "ksmooth", "smoothSpl"),
         nodes=NULL, perDose=TRUE) {
    toType     <- match.arg(toType)
    toDoseUnit <- match.arg(toDoseUnit)
    
    ## copy old DVH structure and convert DVH as well as other dose info
    DVH <- x
    if((toType == "asis") && (toDoseUnit != "asis")) {
        ## just change dose unit in DVH and remaining dose values
        cf <- if( (toupper(x$doseUnit) == "CGY") && (toDoseUnit == "GY")) {
            1/100
        } else if((toupper(x$doseUnit) == "GY")  && (toDoseUnit == "CGY")) {
            100
        } else {
            NA_real_
        }

        DVH$dvh[ , "dose"] <- cf*x$dvh[ , "dose"]
        DVH$doseMin   <- cf*x$doseMin
        DVH$doseMax   <- cf*x$doseMax
        DVH$doseRx    <- cf*x$doseRx
        DVH$isoDoseRx <- cf*x$isoDoseRx
        DVH$doseAvg   <- cf*x$doseAvg
        DVH$doseMed   <- cf*x$doseMed
        DVH$doseMode  <- cf*x$doseMode
        DVH$doseSD    <- cf*x$doseSD
        DVH$doseUnit  <- toDoseUnit
    } else if((toType != "asis") && (toDoseUnit == "asis")) {
        ## just change DVH type
        if(toType == "differential") {
            ## from cumulative to differential
            DVH$dvhDiff <- convertDVH(x$dvh, toType=toType, toDoseUnit="asis",
                                      interp=interp, nodes=nodes, perDose=perDose)
        } else {
            ## from differential to cumulative
            DVH$dvh <- convertDVH(x$dvhDiff, toType=toType, toDoseUnit="asis",
                                  interp=interp, nodes=nodes, perDose=perDose)
        }

        DVH$DVHtype <- toType
    } else if((toType != "asis") && (toDoseUnit != "asis")) {
        ## change DVH type and dose unit
        cf <- if( (toupper(x$doseUnit) == "CGY") && (toDoseUnit == "GY")) {
            1/100
        } else if((toupper(x$doseUnit) == "GY")  && (toDoseUnit == "CGY")) {
            100
        } else {
            NA_real_
        }

        if(toType == "differential") {
            ## from cumulative to differential
            DVH$dvhDiff <- convertDVH(x$dvh, toType=toType, toDoseUnit=toDoseUnit,
                                      interp=interp, nodes=nodes, perDose=perDose)
        } else {
            ## from differential to cumulative
            DVH$dvh <- convertDVH(x$dvhDiff, toType=toType, toDoseUnit=toDoseUnit,
                                  interp=interp, nodes=nodes, perDose=perDose)
        }

        DVH$doseMin   <- cf*x$doseMin
        DVH$doseMax   <- cf*x$doseMax
        DVH$doseRx    <- cf*x$doseRx
        DVH$isoDoseRx <- cf*x$isoDoseRx
        DVH$doseAvg   <- cf*x$doseAvg
        DVH$doseMed   <- cf*x$doseMed
        DVH$doseMode  <- cf*x$doseMode
        DVH$doseSD    <- cf*x$doseSD
        DVH$DVHtype   <- toType
        DVH$doseUnit  <- toDoseUnit
    }
    
    return(DVH)
}

## for list of DVH objects
convertDVH.DVHLst <-
function(x, toType=c("asis", "cumulative", "differential"),
         toDoseUnit=c("asis", "GY", "CGY"),
         interp=c("asis", "linear", "spline", "ksmooth", "smoothSpl"),
         nodes=NULL, perDose=TRUE) {
    toType     <- match.arg(toType)
    toDoseUnit <- match.arg(toDoseUnit)
    interp     <- match.arg(interp)

    dvhL <- Map(convertDVH, x, toType=toType, toDoseUnit=toDoseUnit,
                interp=interp, nodes=list(nodes), perDose=perDose)
    names(dvhL) <- names(x)
    class(dvhL) <- "DVHLst"
    attr(dvhL, which="byPat") <- attributes(x)$byPat

    return(dvhL)
}

## for list of list of DVH objects
convertDVH.DVHLstLst <-
function(x, toType=c("asis", "cumulative", "differential"),
         toDoseUnit=c("asis", "GY", "CGY"),
         interp=c("asis", "linear", "spline", "ksmooth", "smoothSpl"),
         nodes=NULL, perDose=TRUE) {
    toType     <- match.arg(toType)
    toDoseUnit <- match.arg(toDoseUnit)
    interp     <- match.arg(interp)

    dvhLL <- Map(convertDVH, x, toType=toType, toDoseUnit=toDoseUnit,
                 interp=interp, nodes=list(nodes), perDose=perDose)
    names(dvhLL) <- names(x)
    class(dvhLL) <- "DVHLstLst"
    attr(dvhLL, which="byPat") <- attributes(x)$byPat

    return(dvhLL)
}
