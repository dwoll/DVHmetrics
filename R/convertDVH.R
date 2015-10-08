## convert DVH to cumulative/differential DVH
convertDVH <-
function(x, toType=c("asis", "cumulative", "differential"),
         toDoseUnit=c("asis", "GY", "CGY"),
         interp=c("asis", "linear"),#, "spline", "ksmooth", "smoothSpl"),
         nodes=NULL, perDose=TRUE) {
    UseMethod("convertDVH")
}

## for DVH matrix itself
convertDVH.matrix <-
function(x, toType=c("asis", "cumulative", "differential"),
         toDoseUnit=c("asis", "GY", "CGY"),
         interp=c("asis", "linear"),#, "spline", "ksmooth", "smoothSpl"),
         nodes=NULL, perDose=TRUE) {
    toType     <- match.arg(toType)
    toDoseUnit <- match.arg(toDoseUnit)
    interp     <- match.arg(interp)
    
    ## matrix may include duplicated rows
    x <- unique(x)

    ## when dose entries are duplicated with different volumes,
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

    ## number of DVH points
    N <- nrow(x)
    
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
            sm <- switch(interp,
                         ksmooth=getKSmooth(    doseConv, doseRel, volume, volumeRel, nodes=nodes),
                         smoothSpl=getSmoothSpl(doseConv, doseRel, volume, volumeRel, nodes=nodes),
                         getInterpSpl(          doseConv, doseRel, volume, volumeRel, nodes=nodes)) # default

            ## TODO: re-normalize
            doseNew      <- sm[ , "dose"]
            doseRelNew   <- sm[ , "doseRel"]
            volumeNew    <- sm[ , "volume"]
            volumeRelNew <- sm[ , "volumeRel"]
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
            sm <- switch(interp,
                         ksmooth=getKSmooth(    doseConv, doseRel, volume, volumeRel, nodes=nodes),
                         smoothSpl=getSmoothSpl(doseConv, doseRel, volume, volumeRel, nodes=nodes),
                         getInterpSpl(          doseConv, doseRel, volume, volumeRel, nodes=nodes)) # default

            ## TODO: re-normalize
            doseConv  <- sm[ , "dose"]
            doseRel   <- sm[ , "doseRel"]
            volume    <- sm[ , "volume"]
            volumeRel <- sm[ , "volumeRel"]
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
        
        ## differential DVH -> volume is per Gy -> mult with bin-width
        binW <- diff(c(-doseConv[1], doseConv))
        volumeBin <- if(perDose) {
            volume*binW
        } else {
            volume
        }

        volumeRelBin <- if(perDose) {
            volumeRel*binW
        }  else {
            volumeRel
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
            sm <- switch(interp,
                  ksmooth=getKSmooth(  doseNew, doseRelNew, volumeNew, volumeRelNew, nodes=nodes),
                smoothSpl=getSmoothSpl(doseNew, doseRelNew, volumeNew, volumeRelNew, nodes=nodes),
                getInterpSpl(          doseNew, doseRelNew, volumeNew, volumeRelNew, nodes=nodes)) # default

            ## TODO: re-normalize
            doseNew      <- sm[ , "dose"]
            doseRelNew   <- sm[ , "doseRel"]
            volumeNew    <- sm[ , "volume"]
            volumeRelNew <- sm[ , "volumeRel"]
        }
    }

    cbind(dose=doseNew, doseRel=doseRelNew,
          volume=volumeNew, volumeRel=volumeRelNew)
}

## for 1 DVH object
convertDVH.DVHs <-
function(x, toType=c("asis", "cumulative", "differential"),
         toDoseUnit=c("asis", "GY", "CGY"),
         interp=c("asis", "linear"),#, "spline", "ksmooth", "smoothSpl"),
         nodes=NULL, perDose=TRUE) {
    toType     <- match.arg(toType)
    toDoseUnit <- match.arg(toDoseUnit)
    
    ## copy old DVH structure and convert DVH as well as other dose info
    DVH <- x
    if(((toType     == "asis") || (toType     == x$DVHtype)) &&
       ((toDoseUnit != "asis") && (toDoseUnit != x$doseUnit))) {
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
    } else if(((toType     != "asis") && (toType     != x$DVHtype)) &&
              ((toDoseUnit == "asis") || (toDoseUnit == x$doseUnit))) {
        ## just change DVH type
        DVHout <- convertDVH(x$dvh, toType=toType, toDoseUnit="asis",
                             interp=interp, nodes=nodes, perDose=perDose)

        if(toType == "differential") {
            ## if differential -> use dvhDiff and copy cumulative DVH
            DVH$dvh     <- x$dvh
            DVH$dvhDiff <- DVHout
        } else {
            ## if cumulative -> use dvh and copy differential DVH
            DVH$dvh     <- DVHout
            DVH$dvhDiff <- x$dvh
        }

        DVH$DVHtype <- toType
    } else if(((toType     != "asis") && (toType     != x$DVHtype)) &&
              ((toDoseUnit != "asis") && (toDoseUnit != x$doseUnit))) {
        ## change DVH type and dose unit
        cf <- if( (toupper(x$doseUnit) == "CGY") && (toDoseUnit == "GY")) {
            1/100
        } else if((toupper(x$doseUnit) == "GY")  && (toDoseUnit == "CGY")) {
            100
        } else {
            NA_real_
        }

        DVHout <- convertDVH(x$dvh, toType=toType, toDoseUnit=toDoseUnit,
                             interp=interp, nodes=nodes, perDose=perDose)

        if(toType == "differential") {
            ## if differential -> use dvhDiff and copy cumulative DVH
            DVH$dvh     <- x$dvh
            DVH$dvhDiff <- DVHout
        } else {
            ## if cumulative -> use dvh and copy differential DVH
            DVH$dvh     <- DVHout
            DVH$dvhDiff <- x$dvh
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
         interp=c("asis", "linear"),#, "spline", "ksmooth", "smoothSpl"),
         nodes=NULL, perDose=TRUE) {
    toType     <- match.arg(toType)
    toDoseUnit <- match.arg(toDoseUnit)
    interp     <- match.arg(interp)

    dvhL <- Map(convertDVH, x, toType=toType, toDoseUnit=toDoseUnit,
                interp=interp, nodes=list(nodes), perDose=list(perDose))
    names(dvhL) <- names(x)
    class(dvhL) <- "DVHLst"
    attr(dvhL, which="byPat") <- attributes(x)$byPat

    return(dvhL)
}

## for list of list of DVH objects
convertDVH.DVHLstLst <-
function(x, toType=c("asis", "cumulative", "differential"),
         toDoseUnit=c("asis", "GY", "CGY"),
         interp=c("asis", "linear"),#, "spline", "ksmooth", "smoothSpl"),
         nodes=NULL, perDose=TRUE) {
    toType     <- match.arg(toType)
    toDoseUnit <- match.arg(toDoseUnit)
    interp     <- match.arg(interp)

    dvhLL <- Map(convertDVH, x, toType=toType, toDoseUnit=toDoseUnit,
                 interp=interp, nodes=list(nodes), perDose=list(perDose))
    names(dvhLL) <- names(x)
    class(dvhLL) <- "DVHLstLst"
    attr(dvhLL, which="byPat") <- attributes(x)$byPat

    return(dvhLL)
}
