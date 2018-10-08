#####---------------------------------------------------------------------------
## parse character vector from PRIMO DVH file
parsePRIMO <- function(x, planInfo=FALSE, courseAsID=FALSE, ...) {
    planInfo <- as.character(planInfo)

    ## function to extract one information element from a number of lines
    ## make sure only first : is matched -> not greedy
    getElem <- function(pattern, ll, trim=TRUE, iCase=FALSE, collWS=TRUE) {
        line <- ll[grep(pattern, ll)]
        elem <- sub("^.+?:[[:blank:]]+([[:alnum:][:punct:][:blank:]]+$)", "\\1",
                    line, ignore.case=iCase, perl=TRUE)
        elem <- if(trim) {
            trimWS(elem, side="both")
        } else {
            elem
        }

        if(collWS) {
            collWS(elem)
        } else {
            elem
        }
    }

    getDose <- function(pattern, ll, doseRx) {
        line <- trimWS(ll[grep(pattern, ll)], side="both")
        pat  <- "^.+?:[^[:digit:]]*([[:digit:][:punct:]]+) \\([[:alnum:][:punct:]]+\\)$"
        elem <- gsub(pat, "\\1", line, perl=TRUE, ignore.case=TRUE)
        num  <- as.numeric(trimWS(elem))
        if(is.na(num)) {
            warning("No dose found")
        }

        if(grepl("\\[%\\]", line)) {
            ## relative dose
            if(!missing(doseRx)) {
                doseRx * num/100
            } else {
                NA_real_
            }
        } else {
            num
        }
    }

    getDoseUnit <- function(ll) {
        line <- ll[grep("^# Ave\\. Dose.+:", ll)]
        elem <- sub("^.+\\[(GY|CGY|EV/G)\\][[:blank:]]*:.+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        toupper(trimWS(elem))
    }

    getVolUnit <- function(ll) {
        line <- ll[grep("^# Volume.+:", ll)]
        elem <- sub("^.+\\[(CM.+)\\][[:blank:]]*:.+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        toupper(trimWS(elem))
    }

    getDVHtype <- function(ll) {
        line <- ll[grep("^# Type:", ll)]
        elem <- sub("^# Type:[[:blank:]]+(Cumulative|Differential)", "\\1",
                    line, perl=TRUE, ignore.case=TRUE)
        tolower(trimWS(elem))
    }

    getDVHdate <- function(ll) {
        line <- ll[grep("^# Date:", ll)]
        elem <- substr(sub("^# Date: (.+)$", "\\1", line), start=1, stop=10)
        as.Date(elem, format="%d.%m.%Y")
    }

    sStart <- grep("^# Structure: ", x)     # start sections
    sLen   <- diff(c(sStart, length(x)+1))  # length of sections
    if((length(sLen) < 1L) || all(sLen < 1L)) {
        stop("No structures found")
    }

    structList <- Map(function(start, len) x[start:(start+len-1)], sStart, sLen)

    ## extract file header and header info
    header  <- x[seq_len(sStart[1]-2)]                        # header
    patName <- gsub("[^a-z0-9]", "\\1", tempfile(pattern="", tmpdir=""))
    patID   <- gsub("[^a-z0-9]", "\\1", tempfile(pattern="", tmpdir=""))
    plan    <- getElem("^# Project:",  header)
    DVHdate <- getDVHdate(header)
    DVHtype <- getDVHtype(header)

    ## extract DVH from one structure section and store in a list
    ## with DVH itself as a matrix
    getDVH <- function(strct, info) {
        doseRx     <- NA_real_
        doseRxUnit <- NA_character_
        isoDoseRx  <- NA_real_
        doseUnit   <- getDoseUnit(strct)

        ## extract structure, volume, dose min, max, mean, median and sd
        structure <- getElem("^# Structure:", strct)
        structVol <- as.numeric(getElem("^# Volume.+:", strct))
        doseMin   <- getDose("^# Min\\. Dose.+:", strct, doseRx)
        doseMax   <- getDose("^# Max\\. Dose.+:", strct, doseRx)
        doseAvg   <- getDose("^# Ave\\. Dose.+:", strct, doseRx)
        doseMed   <- NA_real_
        doseMode  <- NA_real_
        doseSD    <- NA_real_

        volumeUnit <- getVolUnit(strct)
        volumeUnit <- if(grepl("^CM.+", volumeUnit)) {
            "CC"
        } else if(grepl("^%", volumeUnit)) {
            "PERCENT"
        } else {
            warning("Could not determine volume measurement unit")
            NA_character_
        }

        ## find DVH
        ## DVH column headers
        colHead  <- grep("^# \\.+$", strct)[1] - 1  # column header
        dvhStart <- colHead+2                       # first numeric line of DVH
        dvhEnd   <- grep("^# \\.+$", strct)[2] - 1  # last numeric line of DVH
        dvhLen   <- dvhEnd - dvhStart + 1           # number of DVH lines
        if((length(dvhLen) < 1L) || dvhLen < 1L) {
            stop("No DVH data found")
        }

        ## column headers
        vars1 <- unlist(strsplit(sub("^# (.+)$", "\\1", strct[colHead]),
                                 split="\t"))

        ## remove leading and trailing white space
        vars2 <- tolower(trimWS(vars1))

        ## make sure we recognize all columns in the DVH
        patDose    <- "^dose"
        patDoseRel <- "^rel\\. dose"
        patVol     <- "^volume \\[[^%]\\]"
        patVolRel  <- "^volume \\[%\\]"
        patVolD    <- "^dv/dd"  # differential DVH
        hits <- sum(c(grepl(patDose, vars2), grepl(patDoseRel, vars2),
                      grepl(patVol,  vars2), grepl(patVolRel,  vars2),
                      grepl(patVolD, vars2)))
        if(length(vars2) != hits) {
            stop(c("Could not identify all DVH columns"),
                 paste(vars2, collapse=", "))
        }

        ## replace column headers
        vars3 <- vars2
        vars3[grep(patDose,    vars2)] <- "dose"
        vars3[grep(patDoseRel, vars2)] <- "doseRel"
        vars3[grep(patVol,     vars2)] <- "volume"
        vars3[grep(patVolRel,  vars2)] <- "volumeRel"
        vars3[grep(patVolD,    vars2)] <- "volumeRel"

        ## extract DVH as a matrix and store preceding information
        ## check if dvh is all blank -> no data
        if(all(!nzchar(strct[dvhStart:dvhEnd]))) {
            return(NULL)
        }

        dvh <- data.matrix(read.table(text=strct[dvhStart:dvhEnd],
                                      header=FALSE, stringsAsFactors=FALSE,
                                      colClasses=rep("numeric", length(vars3)),
                                      comment.char="", nrows=dvhLen))

        ## rename
        colnames(dvh) <- vars3

        ## catch special case: structVol is 0.0 due to limited precision
        if("volume" %in% vars3) {
            structVol <- if(info$DVHtype == "cumulative") {
                max(c(structVol, dvh[ , "volume"]))
            } else {
                ## reconstruct volumes
                max(c(structVol, sum(dvh[ , "volume"])))
            }
        }

        ## if we have relative and absolute dose: add doseRx
        if(("doseRel" %in% vars3) && ("dose" %in% vars3)) {
            doseRx <- round(mean(dvh[ , "dose"] * 100 / dvh[ , "doseRel"],
                                 na.rm=TRUE), digits=2)
        }

        ## add information we don't have yet: relative/absolute volume
        if((       "volumeRel" %in% vars3) && !("volume"    %in% vars3)) {
            dvh <- cbind(dvh, volume=structVol*(dvh[ , "volumeRel"]/100))
        } else if(("volume"    %in% vars3) && !("volumeRel" %in% vars3)) {
            dvh <- cbind(dvh, volumeRel=100*(dvh[ , "volume"]/structVol))
        }

        ## add information we don't have yet: relative/absolute dose
        ## considering isoDoseRx
        if((    "doseRel" %in% vars3) && !("dose"    %in% vars3)) {
            if(is.na(isoDoseRx)) {
                warning("I assume that iso-dose is 1")
                dvh <- cbind(dvh, dose=dvh[ , "doseRel"]*doseRx)
            } else {
                dvh <- cbind(dvh, dose=dvh[ , "doseRel"]*doseRx / isoDoseRx)
                # (doseRx/(isoDoseRx/100))*(dvh$doseRel/100)
            }
        } else if(("dose" %in% vars3) && !("doseRel" %in% vars3)) {
            if(is.na(isoDoseRx)) {
                warning("I assume that iso-dose is 1")
                dvh <- cbind(dvh, doseRel=dvh[ , "dose"] / doseRx)
            } else {
                dvh <- cbind(dvh, doseRel=dvh[ , "dose"]*isoDoseRx / doseRx)
                # 100*(dvh$dose/(doseRx/(isoDoseRx/100)))
            }
        }

        ## check if dose is increasing
        stopifnot(isIncreasing(dvh))

        DVH <- list(dvh=dvh,
                    patName=info$patName,
                    patID=info$patID,
                    date=info$date,
                    DVHtype=info$DVHtype,
                    plan=info$plan,
                    course=info$course,
                    quadrant=info$quadrant,
                    structure=structure,
                    structVol=structVol,
                    doseUnit=doseUnit,
                    volumeUnit=volumeUnit,
                    doseMin=doseMin,
                    doseMax=doseMax,
                    doseRx=doseRx,
                    doseRxUnit=doseRxUnit,
                    isoDoseRx=isoDoseRx,
                    doseAvg=doseAvg,
                    doseMed=doseMed,
                    doseMode=doseMode,
                    doseSD=doseSD)

        ## convert differential DVH to cumulative
        ## and add differential DVH separately
        if(info$DVHtype == "differential") {
            DVH$dvh     <- convertDVH(dvh, toType="cumulative",
                                      toDoseUnit="asis", perDose=FALSE)
            DVH$dvhDiff <- dvh
        }

        ## set class
        class(DVH) <- "DVHs"
        return(DVH)
    }

    ## list of DVH data frames with component name = structure
    info <- list(patID=patID, patName=patName, date=DVHdate,
                 plan=plan, DVHtype=DVHtype)
    dvhL <- lapply(structList, getDVH, info=info)
    dvhL <- Filter(Negate(is.null), dvhL)
    names(dvhL) <- sapply(dvhL, function(y) y$structure)
    if(length(unique(names(dvhL))) < length(dvhL)) {
        warning("Some structures have the same name - this can lead to problems")
    }

    class(dvhL) <- "DVHLst"
    attr(dvhL, which="byPat") <- TRUE

    return(dvhL)
}
