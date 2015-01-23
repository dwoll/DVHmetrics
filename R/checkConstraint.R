## S3 generic function
checkConstraint <-
function(x, constr, byPat=TRUE, semSign=FALSE,
         sortBy=c("none", "observed", "compliance", "structure",
                  "constraint", "patID", "deltaV", "deltaD", "dstMin")) {
    UseMethod("checkConstraint")
}

checkConstraint.DVHs <-
function(x, constr, byPat=TRUE, semSign=FALSE,
         sortBy=c("none", "observed", "compliance", "structure",
                  "constraint", "patID", "deltaV", "deltaD", "dstMin")) {
    x <- if(byPat) {
        setNames(list(x), x$patID)
    } else {
        setNames(list(x), x$structure)
    }

    class(x) <- c("DVHLst", "list")
    attr(dvhL, which="byPat") <- byPat

    #NextMethod("checkConstraint")
    checkConstraint.DVHLst(x, constr=constr, byPat=byPat, semSign=semSign,
                            sortBy=sortBy)
}

## with byPat=TRUE
## x is DVH list for 1 id,        constr is list of constraints, 1 per structure
## with byPat=FALSE
## x is DVH list for 1 structure, constr is list of constraints, 1 per id
checkConstraint.DVHLst <-
function(x, constr, byPat=TRUE, semSign=FALSE,
         sortBy=c("none", "observed", "compliance", "structure",
                  "constraint", "patID", "deltaV", "deltaD", "dstMin")) {
    ## make sure DVH list is organized as required for byPat
    if(is.null(attributes(x)$byPat) || attributes(x)$byPat != byPat) {
        stop(c("DVH list organization by-patient / by-structure",
               " either could not be determined or is different from byPat"))
    }

    xConstrSub  <- harmoConstrDVH(x, constr=constr, byPat=byPat)
    constrParse <- Map(parseConstraint, xConstrSub$constr)

    ## calculate difference between 1 observed metric value and 1 constraint
    cmpMetrics <- function(observed, smInv, DV, valCmp, cmp, valRef, dstInfo) {
        valDiff     <- observed - valCmp
        valRelPC    <- 100 * (valDiff/valCmp)
        valDiffInv  <- smInv - valRef
        valRelPCInv <- 100 * (valDiffInv/valRef)

        deltaV   <- if(DV == "V") { valDiff  } else if(DV == "D") { valDiffInv  }
        deltaVpc <- if(DV == "V") { valRelPC } else if(DV == "D") { valRelPCInv }
        deltaD   <- if(DV == "D") { valDiff  } else if(DV == "V") { valDiffInv  }
        deltaDpc <- if(DV == "D") { valRelPC } else if(DV == "V") { valRelPCInv }

        list(observed=observed, deltaV=deltaV, deltaVpc=deltaVpc,
             deltaD=deltaD, deltaDpc=deltaDpc, valCmp=valCmp,
             cmp=cmp, dstMin=dstInfo$dstMin, ptMinD=dstInfo$ptMinD,
             ptMinV=dstInfo$ptMinV)
    }

    ## calculate metric and corresponding inverse metrics from given
    ## constraint set - vectorized in metrics
    stridMetrics <- function(dvh, cnstr) {
        observed <- ifelse(cnstr$valid,
                           unlist(getMetric(dvh, metric=cnstr$metric)),
                           NA_real_)

        
        ## inverse metrics
        smInv <- ifelse(cnstr$valid,
                        unlist(getMetric(dvh, metric=cnstr$metricInv)),
                        NA_real_)

        names(observed) <- cnstr$constraint
        names(smInv)    <- cnstr$constraint

        ## distance constraint to closest point on DVH curve
        Dcoord  <- with(cnstr, ifelse(DV == "D", valCmp, valRef))
        Vcoord  <- with(cnstr, ifelse(DV == "V", valCmp, valRef))
        Dcoord  <- ifelse(cnstr$valid, Dcoord, NA_real_)
        Vcoord  <- ifelse(cnstr$valid, Vcoord, NA_real_)
        volRel  <- with(cnstr, ((DV == "D") & (unitRef == "%")) |
                               ((DV == "V") & (unitCmp == "%")))
        doseRel <- with(cnstr, ((DV == "D") & (unitCmp == "%")) |
                               ((DV == "V") & (unitRef == "%")))

        dstDVH <- dvhDistance(dvh,
                              DV=data.frame(D=Dcoord, V=Vcoord,
                                            volRel=volRel, doseRel=doseRel,
                                            unitRef=cnstr$unitRef, unitCmp=cnstr$unitCmp))

        with(cnstr, Map(cmpMetrics, observed=observed, smInv=smInv,
                        DV=DV, valCmp=valCmp, cmp=cmp, valRef=valRef,
                        dstInfo=dstDVH))
    }

    ## get requested metrics for each structure/id
    metrics <- Map(stridMetrics, xConstrSub$x, constrParse)

    ## check constraints
    leqGeq <- function(y) {
        observed <- y$observed
        valCmp   <- y$valCmp
        fun      <- y$cmp
        do.call(fun, list(observed, valCmp))
    }

    cmpFun <- function(y) {
        Map(leqGeq, y)
    }

    compL <- Map(cmpFun, metrics)
    metrL <- melt(metrics, value.name="observed")

    ## transform into data frame
    compDF <- melt(compL, value.name="compliance")
    metrDF <- dcast(metrL, L1 + L2 ~ L3, value.var="observed")
    resDF  <- merge(compDF, metrDF)
    names(resDF)[names(resDF) == "L1"] <- if(byPat) {
        "structure"
    } else {
        "patID"
    }

    names(resDF)[names(resDF) == "L2"] <- "constraint"

    ## make sign of deltaD/deltaV semantically indicate compliance
    ## negative -> compliance, positive -> no compliance
    if(semSign) {
        resDF$deltaD   <- ifelse(resDF$compliance,
                                 -1*resDF$deltaD  *sign(resDF$deltaD),
                                    resDF$deltaD  *sign(resDF$deltaD))
        resDF$deltaDpc <- ifelse(resDF$compliance,
                                 -1*resDF$deltaDpc*sign(resDF$deltaDpc),
                                    resDF$deltaDpc*sign(resDF$deltaDpc))
        resDF$deltaV   <- ifelse(resDF$compliance,
                                 -1*resDF$deltaV  *sign(resDF$deltaV),
                                    resDF$deltaV  *sign(resDF$deltaV))
        resDF$deltaVpc <- ifelse(resDF$compliance,
                                 -1*resDF$deltaVpc*sign(resDF$deltaVpc),
                                    resDF$deltaVpc*sign(resDF$deltaVpc))
        resDF$dstMin   <- ifelse(resDF$compliance,
                                 -1*resDF$dstMin  *sign(resDF$dstMin),
                                    resDF$dstMin  *sign(resDF$dstMin))
    }

    ## reorder variables to return
    if(byPat) {
        patID     <- xConstrSub$x[[1]]$patID
    } else {
        structure <- xConstrSub$x[[1]]$structure
    }
        
    finDF <- with(resDF,
                  data.frame(patID, structure, constraint, observed, compliance,
                             deltaV, deltaVpc, deltaD, deltaDpc, dstMin, ptMinD, ptMinV,
                             stringsAsFactors=FALSE))

    ## sort
    if(!("none" %in% sortBy)) {
        sortBySub <- sortBy[sortBy %in% names(finDF)]
        idx       <- do.call("order", lapply(sortBySub, function(y) with(finDF, get(y))))
        finDF     <- finDF[idx, ]
    }

    return(finDF)
}

## x is a DVH list (1 per id/structure) of lists
## constr is a list (1 per id/structure) of lists (1 for each id/structure)
checkConstraint.DVHLstLst <-
function(x, constr, byPat=TRUE, semSign=FALSE,
         sortBy=c("none", "observed", "compliance", "structure",
                  "constraint", "patID", "deltaV", "deltaD", "dstMin")) {

    ## re-organize x into by-patient or by-structure form if necessary
    isByPat <- attributes(x)$byPat
    xRO <- if(is.null(isByPat) || (isByPat != byPat)) {
        reorgByPat(x, byPat=byPat)
    } else {
        x
    }

    ## make sure DVH and constraint have the same IDs / structures
    xConstrSub <- harmoConstrDVH(xRO, constr=constr, byPat=byPat)

    ## check the constraints
    res <- Map(checkConstraint,
               x=xConstrSub$x, constr=xConstrSub$constr,
               byPat=byPat, semSign=semSign)

    ## transform into data frame
    resL    <- melt(res, id.vars=c("patID", "structure", "constraint", "compliance"))
    resL$L1 <- NULL   # redundant with patID (bPat=TRUE) or structure (byPat=FALSE)
    resDF   <- dcast(resL, patID + structure + constraint + compliance ~ variable,
                     value.var="value")

    ## reorder variables to return
    finDF <- with(resDF, data.frame(patID, structure, constraint, observed, compliance,
                                    deltaV, deltaVpc, deltaD, deltaDpc,
                                    dstMin, ptMinD, ptMinV, stringsAsFactors=FALSE))

    ## sort
    if(!("none" %in% sortBy)) {
        sortBySub <- sortBy[sortBy %in% names(finDF)]
        idx <- do.call("order", lapply(sortBySub, function(y) with(finDF, get(y))))
        finDF <- finDF[idx, ]
    }

    return(finDF)
}
