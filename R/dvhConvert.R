## convert DVH to cumulative/differential DVH
dvhConvert <-
function(x, toType=c("asis", "cumulative", "differential"),
         toDoseUnit=c("asis", "GY", "CGY")) {
    UseMethod("dvhConvert")
}

## for DVH matrix itself
dvhConvert.matrix <-
function(x, toType=c("asis", "cumulative", "differential"),
         toDoseUnit=c("asis", "GY", "CGY")) {
    toType     <- match.arg(toType)
    toDoseUnit <- match.arg(toDoseUnit)

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
        doseNew      <- doseConv
        doseRelNew   <- doseRel
        volumeNew    <- volume
        volumeRelNew <- volumeRel
    } else if(toType == "cumulative") {
        ## convert from differential to cumulative DVH
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
        volumeBin    <- volume   *diff(c(-doseConv[1], doseConv))
        volumeRelBin <- volumeRel*diff(c(-doseConv[1], doseConv))

        ## volume "survival" curve starting at max volume
        volumeNew    <- sum(volumeBin)    - c(0, cumsum(volumeBin))
        volumeNewRel <- sum(volumeRelBin) - c(0, cumsum(volumeRelBin))
    } else if(toType == "differential") {
        ## convert from cumulative to differential DVH
        ## dose category half-widths
        doseCatHW    <- diff(doseConv)/2
        doseRelCatHW <- diff(doseRel)/2

        doseNew      <- doseConv[-length(doseConv)] + doseCatHW
        doseRelNew   <-  doseRel[-length(doseRel)]  + doseRelCatHW

        ## differential DVH -> volume is per Gy -> divide by bin-width
        volumeNew    <- -diff(volume)    / (2*doseCatHW)
        volumeNewRel <- -diff(volumeRel) / (2*doseRelCatHW)
    }

    cbind(dose=doseNew, doseRel=doseRelNew,
          volume=volumeNew, volumeRel=volumeNewRel)
}

## for 1 DVH object
dvhConvert.DVHs <-
function(x, toType=c("asis", "cumulative", "differential"),
         toDoseUnit=c("asis", "GY", "CGY")) {
    toType     <- match.arg(toType)
    toDoseUnit <- toupper(match.arg(toDoseUnit))
    
    ## copy old DVH structure and convert DVH as well as other dose info
    DVH <- x
    if(((toType     == "asis") || (toType     == x$DVHtype)) &&
       ((toDoseUnit != "ASIS") && (toDoseUnit != x$doseUnit))) {
        ## just change dose unit in DVH and remaining dose values
        cf <- if(toDoseUnit == "GY") { 1/100 } else if(toDoseUnit == "CGY") { 100 }
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
              ((toDoseUnit == "ASIS") || (toDoseUnit == x$doseUnit))) {
        ## just change DVH type
        DVH$dvh  <- dvhConvert(x$dvh, toType=toType, toDoseUnit="asis")
        DVH$DVHtype <- toType
    } else if(((toType     != "asis") && (toType     != x$DVHtype)) &&
              ((toDoseUnit != "ASIS") && (toDoseUnit != x$doseUnit))) {
        ## change DVH type and dose unit
        cf <- if(toDoseUnit == "GY") { 1/100 } else if(toDoseUnit == "CGY") { 100 }
        DVH$dvh <- dvhConvert(x$dvh, toType=toType, toDoseUnit=toDoseUnit)
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
dvhConvert.DVHLst <-
function(x, toType=c("asis", "cumulative", "differential"),
         toDoseUnit=c("asis", "GY", "CGY")) {
    toType     <- match.arg(toType)
    toDoseUnit <- match.arg(toDoseUnit)

    dvhL <- Map(dvhConvert, x, toType=toType, toDoseUnit=toDoseUnit)
    names(dvhL) <- names(x)
    class(dvhL) <- c("DVHLst", "list")
    attr(dvhL, which="byPat") <- attributes(x)$byPat

    return(dvhL)
}

## for list of list of DVH objects
dvhConvert.DVHLstLst <-
function(x, toType=c("asis", "cumulative", "differential"),
         toDoseUnit=c("asis", "GY", "CGY")) {
    toType     <- match.arg(toType)
    toDoseUnit <- match.arg(toDoseUnit)

    dvhLL <- Map(dvhConvert, x, toType=toType, toDoseUnit=toDoseUnit)
    names(dvhLL) <- names(x)
    class(dvhLL) <- c("DVHLstLst", "list")
    attr(dvhLL, which="byPat") <- attributes(x)$byPat

    return(dvhLL)
}
