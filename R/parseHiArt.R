#####---------------------------------------------------------------------------
## parse character vector from Tomo HiArt DVH file
parseHiArt <- function(x, planInfo=FALSE, courseAsID=FALSE, ...) {
    planInfo <- as.character(planInfo)

    ## find columns for structure, dose, volume
    vars <- as.matrix(read.csv(text=x[1], header=FALSE,
                               stringsAsFactors=FALSE, comment.char="")[1, ])

    structIdx <- seq(1, length(vars), by=3)
    doseIdx   <- seq(2, length(vars), by=3)
    volumeIdx <- seq(3, length(vars), by=3)

    ## get dose and volume units
    varDose  <- vars[grepl("Dose",   vars)][1]
    varVol   <- vars[grepl("Volume", vars)][1]
    if(grepl("Relative", varDose, ignore.case=TRUE)) {
        ## TODO: need example file for this
        warning("HiArt files with relative dose are not implemented")
        isDoseRel <- TRUE
        doseUnit  <- "PERCENT"
    } else {
        isDoseRel <- FALSE
        doseUnit  <- toupper(sub("Dose \\((.+)\\)", "\\1", varDose))
    }

    if(grepl("Relative", varVol, ignore.case=TRUE)) {
        isVolRel   <- TRUE
        volumeUnit <- "PERCENT"
    } else {
        warning("HiArt files with absolute volume are not implemented")
        ## TODO: need example file for this
        isVolRel   <- FALSE
        volumeUnit <- NA_character_
    }

    ## Tomo HiArt has no patient name or id in file
    ## check if additional information is given via option hiart
    dots <- list(...)
    if(hasName(dots, "hiart")) {
        info <- dots[["hiart"]]
        patName <- if(hasName(info, "patName")) {
            info[["patName"]]
        } else {
            gsub("[^a-z0-9]", "\\1", tempfile(pattern="", tmpdir=""))
        }

        patID <- if(hasName(info, "patID")) {
            info[["patID"]]
        } else {
            gsub("[^a-z0-9]", "\\1", tempfile(pattern="", tmpdir=""))
        }

        doseRx <- if(hasName(info, "doseRx")) {
            info[["doseRx"]]
        } else {
            NA_real_
        }

        structVol <- if(hasName(info, "structVol")) {
            vols <- info[["structVol"]]
            ## expand to all structures
            structures <- sub("(.+)\\(STANDARD\\)$", "\\1", vars[structIdx])
            structVols <- setNames(rep(NA_real_, length(structures)),
                                   structures)
            structVols[names(vols)] <- unlist(vols)
            structVols
        } else {
            NULL
        }

        if(hasName(info, "volumeUnit")) {
            volumeUnit <- info[["volumeUnit"]]
        }
    } else {
        patName   <- gsub("[^a-z0-9]", "\\1", tempfile(pattern="", tmpdir=""))
        patID     <- gsub("[^a-z0-9]", "\\1", tempfile(pattern="", tmpdir="")) # as.character(trunc(runif(1, 0, 10000001)))
        doseRx    <- NA_real_
        structVol <- NA_real_
    }

    plan      <- NA_character_
    isoDoseRx <- NA_real_
    DVHdate   <- NA_character_

    ## read all data
    ## remove all non numbers / delimiters first
    x[-1] <- gsub("[^[:digit:],.]", "", x[-1])
    datAll <- data.matrix(read.csv(text=x[-1], header=FALSE,
                                   stringsAsFactors=FALSE, comment.char=""))

    ## extract DVH from one structure section and store in a list
    ## with DVH itself as a matrix
    getDVH <- function(strIdx, dIdx, vIdx, info, structVol, doseRx) {
        structure <- sub("(.+)\\(STANDARD\\)$", "\\1", vars[strIdx])

        ## extract DVH as a matrix and set variable names
        dvh <- datAll[ , c(dIdx, vIdx)]
        haveVars <- if(isVolRel) {
            if(isDoseRel) {
                c("doseRel", "volumeRel")
            } else {
                c("dose",    "volumeRel")
            }
        } else {
            if(isDoseRel) {
                c("doseRel", "volume")
            } else {
                c("dose",    "volume")
            }
        }

        colnames(dvh) <- haveVars

        ## add information we don't have yet
        ## relative/absolute volume/dose
        if((       "volumeRel" %in% haveVars) && !("volume"    %in% haveVars)) {
            dvh <- cbind(dvh, volume=structVol*(dvh[ , "volumeRel"]/100))
        } else if(("volume"    %in% haveVars) && !("volumeRel" %in% haveVars)) {
            dvh <- cbind(dvh, volumeRel=100*(dvh[ , "volume"]/structVol))
        }

        ## add information we don't have yet: relative/absolute dose
        ## without considering isoDoseRx
        isoDoseRxTmp <- 100
        if((    "doseRel" %in% haveVars) && !("dose"    %in% haveVars)) {
            dvh <- cbind(dvh, dose=dvh[ , "doseRel"]*doseRx / isoDoseRxTmp)
            # (doseRx/(isoDoseRxTmp/100))*(dvh$doseRel/100)
        } else if(("dose" %in% haveVars) && !("doseRel" %in% haveVars)) {
            dvh <- cbind(dvh, doseRel=dvh[ , "dose"]*isoDoseRxTmp / doseRx)
            # 100*(dvh$dose/(doseRx/(isoDoseRxTmp/100)))
        }

        ## check if dose is increasing
        stopifnot(isIncreasing(dvh))

        ## differential or cumulative DVH
        DVHtype <- dvhType(dvh)

        DVH <- list(dvh=dvh,
                    patName=info$patName,
                    patID=info$patID,
                    date=info$date,
                    DVHtype=DVHtype,
                    plan=info$plan,
                    structure=structure,
                    structVol=structVol,
                    doseUnit=info$doseUnit,
                    volumeUnit=info$volumeUnit,
                    doseRx=doseRx,
                    isoDoseRx=isoDoseRx,
                    doseMin=NA_real_,
                    doseMax=NA_real_,
                    doseAvg=NA_real_,
                    doseMed=NA_real_,
                    doseMode=NA_real_,
                    doseSD=NA_real_)

        ## convert differential DVH to cumulative
        ## and add differential DVH separately
        if(DVHtype == "differential") {
            warning("I assume differential DVH is per unit dose\nbut I have no information on this")
            DVH$dvh <- convertDVH(dvh, toType="cumulative",
                                  toDoseUnit="asis", perDose=TRUE)
            DVH$dvhDiff <- dvh
        }

        ## set class
        class(DVH) <- "DVHs"
        return(DVH)
    }

    ## list of DVH data frames with component name = structure
    info <- list(patID=patID, patName=patName, date=DVHdate,
                 plan=plan, doseRx=doseRx, isoDoseRx=isoDoseRx,
                 doseUnit=doseUnit, volumeUnit=volumeUnit)
    dvhL <- Map(getDVH, structIdx, doseIdx, volumeIdx, info=list(info),
                structVol, doseRx)
    dvhL <- Filter(Negate(is.null), dvhL)
    names(dvhL) <- sapply(dvhL, function(y) y$structure)
    if(length(unique(names(dvhL))) < length(dvhL)) {
        warning("Some structures have the same name - this can lead to problems")
    }

    class(dvhL) <- "DVHLst"
    attr(dvhL, which="byPat") <- TRUE

    return(dvhL)
}
