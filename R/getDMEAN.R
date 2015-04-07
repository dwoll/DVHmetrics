getDMEAN <-
function(x, interp=c("linear", "spline", "ksmooth", "smoothSpl")) {
    UseMethod("getDMEAN")
}

getDMEAN.DVHs <-
function(x, interp=c("linear", "spline", "ksmooth", "smoothSpl")) {
    interp <- match.arg(interp)

    ## differential DVH
    xD <- convertDVH(x, toType="differential", interp=interp, nodes=5001L)

    doseMed <- if(interp == "linear") {
        tryCatch(approx(x$dvh[ , "volumeRel"], x$dvh[ , "dose"], 50, method="linear", rule=1)$y,
                 error=function(e) return(NA_real_))
    } else if(interp == "spline") {
        sfun <- try(splinefun(x$dvh[ , "volumeRel"], x$dvh[ , "dose"], method="monoH.FC"))
        if(!inherits(sfun, "try-error")) {
            sfun(50)
        } else {
            NA_real_
        }
    } else if(interp == "ksmooth") {
        bw <- try(KernSmooth::dpill(x$dvh[ , "volumeRel"], x$dvh[ , "dose"]))
        bw <- if(!inherits(bw, "try-error")) {
            bw
        } else {
            NA_real_
        }
        sm <- try(KernSmooth::locpoly(x$dvh[ , "volumeRel"], x$dvh[ , "dose"], bandwidth=bw,
                                  gridsize=5001L, degree=3))
        if(!inherits(sm, "try-error")) {
            idx <- which.min(abs(sm$x-50))
            sm$y[idx]
        } else {
            NA_real_
        }
    } else if(interp == "smoothSpl") {
        sm <- try(smooth.spline(x$dvh[ , "volumeRel"], x$dvh[ , "dose"]))
        if(!inherits(sm, "try-error")) {
            predict(sm, 50)$y
        } else {
            NA_real_
        }
    } else {
        NA_real_
    }

    ## dose category mid-points
    doseMidPt <- xD$dvh[ , "dose"]
    ## differential DVH -> volume is per Gy -> mult with bin-width
    binW      <- diff(c(-xD$dvh[1, "dose"], xD$dvh[ , "dose"]))
    volRelBin <- xD$dvh[ , "volumeRel"]*binW

    ## available volume
    volume <- if(!all(is.na(xD$dvh[ , "volume"]))) {
        xD$dvh[ , "volume"]
    } else {
        xD$dvh[ , "volumeRel"]
    }

    doseMin <- min(xD$dvh[volume > 0, "dose"])
    doseMax <- max(xD$dvh[volume > 0, "dose"])
    doseAvg <- sum(doseMidPt*volRelBin/100)
    ## doseAvg <- sum(doseMidPt*(-1)*diff(x$dvh[, "volumeRel"]/100))
    doseSD  <- sqrt(sum(doseMidPt^2*volRelBin/100) - doseAvg^2)

    ## for mode, abs or rel volume does not matter
    doseMode <- if(!all(is.na(xD$dvh[ , "volume"]))) {
        xD$dvh[which.max(xD$dvh[ , "volume"]),    "dose"]
    } else if(!all(is.na(xD$dvh[ , "volumeRel"]))) {
        xD$dvh[which.max(xD$dvh[ , "volumeRel"]), "dose"]
    } else {
        NA_real_
    }

    doseAvgTPS <- if(!is.null(x$doseAvg)) {
        x$doseAvg
    } else {
        NA_real_
    }

    doseMedTPS <- if(!is.null(x$doseMed)) {
        x$doseMed
    } else {
        NA_real_
    }

    doseMinTPS <- if(!is.null(x$doseMin)) {
        x$doseMin
    } else {
        NA_real_
    }

    doseMaxTPS <- if(!is.null(x$doseMax)) {
        x$doseMax
    } else {
        NA_real_
    }

    metrics <- data.frame(patID=x$patID, structure=x$structure,
                          doseMin=doseMin, doseMax=doseMax, doseAvg=doseAvg,
                          doseMed=doseMed, doseSD=doseSD, doseMode=doseMode,
                          doseAvgTPS=doseAvgTPS, doseMedTPS=doseMedTPS,
                          doseMinTPS=doseMinTPS, doseMaxTPS=doseMaxTPS)

    return(metrics)
}

getDMEAN.DVHLst <-
function(x, interp=c("linear", "spline", "ksmooth", "smoothSpl")) {
    ml <- Map(getDMEAN, x, interp=list(interp))
    df <- do.call("rbind", ml)
    rownames(df) <- NULL
    df
}

getDMEAN.DVHLstLst <-
function(x, interp=c("linear", "spline", "ksmooth", "smoothSpl")) {
    ml <- Map(getDMEAN, x, interp=list(interp))
    df <- do.call("rbind", ml)
    rownames(df) <- NULL
    df
}

## mean from differential DVH from
## integration of dose * (spline fit derivative)
#     xD <- convertDVH(x, toType="differential", smooth=interp, nodes=1001L)
#     ## linear
#     lfun <- approxfun(dose, 100-volume, method="linear", rule=2)
#     meanLin <- tryCatch(pracma::quadgk(function(y) {
#         y*numDeriv::grad(lfun, x=y)/100 }, 0, max(dose)),
#         error=function(e) return(NA_real_))
#
#     ## monotone Hermite spline
#     sfun <- try(splinefun(dose, 100-volume, method="monoH.FC"))
#     meanMHSpl <- if(!inherits(sfun, "try-error")) {
#         tryCatch(pracma::quadgk(function(y) {
#             y*sfun(y, deriv=1)/100 }, 0, max(dose)),
#             error=function(e) return(NA_real_))
#     } else {
#         NA_real_
#     }
#
#     ## locpoly
#     bw <- KernSmooth::dpill(dose, volume)
#     lp <- KernSmooth::locpoly(dose, 100-volume, drv=0, gridsize=10001L,
#                               bandwidth=bw, degree=3)
#     lpfun <- approxfun(lp$x, lp$y, method="linear", rule=2)
#     meanLP <- tryCatch(pracma::quadgk(function(y) {
#         y*numDeriv::grad(lpfun, x=y)/100 }, 0, max(dose)),
#         error=function(e) return(NA_real_))
