#####---------------------------------------------------------------------------
## parse character vector from ProSoma DVH file
parseProSoma <- function(x, planInfo=FALSE, courseAsID=FALSE, ...) {
    planInfo <- as.character(planInfo)

    ## function to extract patient name from first line
    getPatName <- function(ll) {
        paste(trimWS(rev(strsplit(ll, ",")[[1]])), collapse=" ")
    }

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
        line <- ll[grep(pattern, ll)]
        pat  <- "^.+?:[^[:digit:]]*([[:digit:][:punct:]]+)(Gy|cGy|%)*([[:alnum:][:punct:][:blank:]]*)$"
        elem <- gsub(pat, "\\1", line, perl=TRUE, ignore.case=TRUE)
        grp2 <- gsub(pat, "\\2", line, perl=TRUE, ignore.case=TRUE)
        grp3 <- gsub(pat, "\\3", line, perl=TRUE, ignore.case=TRUE)
        if(nzchar(grp2) || nzchar(grp3)) {
            warning("Non-standard dose line found")
        }

        num <- as.numeric(trimWS(elem))
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
        line <- ll[grep("^Average \\[.+\\]:.+", ll)]
        elem <- sub("^.+\\[(GY|CGY)\\]:.+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        toupper(trimWS(elem))
    }

    getVolUnit <- function(ll) {
        line <- ll[grep("^Volume \\[.+\\]:.+", ll)]
        elem <- sub("^.+\\[(.+)\\]:.+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        toupper(trimWS(elem))
    }

    getDVHtype <- function(ll) {
        line <- ll[grep("^(Cumulative|Differential) DVH", ll)]
        elem <- sub("^(Cumulative|Differential).+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        tolower(trimWS(elem))
    }

    sStart <- grep("^Name:[[:blank:]]+[[:alnum:][:punct:]]", x) # start of sections
    sLen   <- diff(c(sStart, length(x)+1))                # length of sections
    if((length(sLen) < 1L) || all(sLen < 1L)) {
        stop("No structures found")
    }

    structList <- Map(function(start, len) x[start:(start+len-1)], sStart, sLen)

    ## extract file header and header info
    ## ProSoma file does not have a PatID
    header     <- x[seq_len(sStart[1]-1)]                        # header
    patName    <- getPatName(header[1])   # patient name
    patID      <- gsub("[^a-z0-9]", "\\1", tempfile(pattern="", tmpdir="")) # as.character(trunc(runif(1, 0, 10000001)))
    plan       <- NA_character_
    doseRx     <- NA_real_
    isoDoseRx  <- NA_real_
    DVHdate    <- NA_character_
    DVHtype    <- getDVHtype(header)
    doseRxUnit <- NA_character_
    doseRx     <- NA_real_

    ## extract DVH from one structure section and store in a list
    ## with DVH itself as a matrix
    getDVH <- function(strct, info) {
        ## extract information from info list
        doseRx     <- info$doseRx
        doseRxUnit <- info$doseRxUnit
        isoDoseRx  <- info$isoDoseRx

        doseUnit <- getDoseUnit(strct)
        if(!grepl("^(GY|CGY)$", doseUnit)) {
            warning("Could not determine dose measurement unit")
            doseUnit <- NA_character_
        }

        ## check if we have dose Rx
        ## if so, does it have the same unit as doseUnit -> convert
        if(!is.na(doseUnit) && !is.na(doseRxUnit)) {
            if((doseUnit == "GY") && (doseRxUnit == "CGY")) {
                doseRx <- doseRx/100
            } else if((doseUnit == "CGY") && (doseRxUnit == "GY")) {
                doseRx <- doseRx*100
            }
        }

        ## extract structure, volume, dose min, max, mean, median and sd
        structure <- getElem("^Name:", strct)
        structVol <- as.numeric(getElem("^Volume \\[.+\\]:.+", strct))
        doseMin   <- getDose("^Min \\[.+\\]:.+",     strct, doseRx)
        doseMax   <- getDose("^Max \\[.+\\]:.+",     strct, doseRx)
        doseAvg   <- getDose("^Average \\[.+\\]:.+", strct, doseRx)
        doseMed   <- getDose("^D50% \\[.+\\]:.+",    strct, doseRx)
        doseMode  <- NA_real_
        doseSD    <- NA_real_

        volumeUnit <- getVolUnit(strct)
        volumeUnit <- if(grepl("^CM.+", volumeUnit)) {
            "CC"
        } else if(grepl("^ML", volumeUnit)) {
            ## ml = cm^3
            "CC"
        } else if(grepl("^%", volumeUnit)) {
            "PERCENT"
        } else {
            warning("Could not determine volume measurement unit")
            NA_character_
        }

        ## find DVH
        ## DVH column headers
        colHead  <- grep("DOSE \\[(%|GY|CGY)\\].+VOLUME", strct,
                         ignore.case=TRUE, perl=TRUE)
        dvhStart <- colHead+1                 # first numeric line of DVH
        dvhLen   <- length(strct) - dvhStart + 1
        if((length(dvhLen) < 1L) || dvhLen < 1L) {
            stop("No DVH data found")
        }

        ## column headers
        vars1 <- unlist(strsplit(strct[colHead],
                        split="\\t+", fixed=FALSE, perl=TRUE))

        ## remove leading and trailing white space
        vars2 <- tolower(trimWS(vars1))

        ## make sure we recognize all columns in the DVH
        patDose    <- "^dose \\[(GY|CGY)\\]"
        patDoseRel <- "^dose \\[%\\]"
        patVolRel  <- "^volume"
        patVolD    <- "^volume"  # differential DVH
        hits <- sum(c(grepl(patDose,    vars2, ignore.case=TRUE),
                      grepl(patDoseRel, vars2, ignore.case=TRUE),
                      grepl(patVolRel,  vars2, ignore.case=TRUE),
                      grepl(patVolD,    vars2, ignore.case=TRUE) &&
                      DVHtype == "differential"))
        if(length(vars2) != hits) {
            stop(c("Could not identify all DVH columns"),
        		 paste(vars2, collapse=", "))
        }

        ## replace column headers
        vars3 <- vars2
        vars3[grep(patDose,    vars2, ignore.case=TRUE)] <- "dose"
        vars3[grep(patDoseRel, vars2, ignore.case=TRUE)] <- "doseRel"
        vars3[grep(patVolRel,  vars2, ignore.case=TRUE)] <- "volumeRel"
        vars3[grepl(patVolD,   vars2, ignore.case=TRUE) &
              (DVHtype == "differential")] <- "volumeRel"

        ## extract DVH as a matrix and store preceding information
        ## read line length(strct) for cases where file does not end with a
        ## blank line -> this will then be last DVH line, otherwise blank
        ## check if dvh is all blank -> no data
        if(all(!nzchar(strct[dvhStart:length(strct)]))) {
            return(NULL)
        }

        dvh <- data.matrix(read.table(text=strct[dvhStart:length(strct)],
                                      header=FALSE, stringsAsFactors=FALSE,
                                      colClasses=rep("numeric", length(vars3)),
                                      comment.char="", nrows=dvhLen))

        ## rename
        colnames(dvh) <- vars3

        ## add information we don't have yet: relative/absolute volume
        if((       "volumeRel" %in% vars3) && !("volume"    %in% vars3)) {
            dvh <- cbind(dvh, volume=structVol*(dvh[ , "volumeRel"]/100))
        } else if(("volume"    %in% vars3) && !("volumeRel" %in% vars3)) {
            dvh <- cbind(dvh, volumeRel=100*(dvh[ , "volume"]/structVol))
        }

        ## add information we don't have yet: relative/absolute dose
        ## considering isoDoseRx
        if((    "doseRel" %in% vars3) && !("dose"    %in% vars3)) {
            dvh <- cbind(dvh, dose=dvh[ , "doseRel"]*doseRx / isoDoseRx)
            # (doseRx/(isoDoseRx/100))*(dvh$doseRel/100)
        } else if(("dose" %in% vars3) && !("doseRel" %in% vars3)) {
            dvh <- cbind(dvh, doseRel=dvh[ , "dose"]*isoDoseRx / doseRx)
            # 100*(dvh$dose/(doseRx/(isoDoseRx/100)))
        }

        ## check if dose is increasing
        stopifnot(isIncreasing(dvh))

        DVH <- list(dvh=dvh,
                    patName=info$patName,
                    patID=info$patID,
                    date=info$date,
                    DVHtype=DVHtype,
                    plan=info$plan,
                    course=NA_character_,
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

        ## convert differential DVH (per unit dose) to cumulative
        ## and add differential DVH separately
        if(info$DVHtype == "differential") {
            # stop("differential DVH currently not supported")
            DVH$dvh     <- convertDVH(dvh,
                                      toType="cumulative",
                                      toDoseUnit="asis",
                                      perDose=FALSE)
            DVH$dvhDiff <- dvh
        }

        ## set class
        class(DVH) <- "DVHs"
        return(DVH)
    }

    ## list of DVH data frames with component name = structure
    info <- list(patID=patID, patName=patName, date=DVHdate,
                 DVHtype=DVHtype, plan=plan, course=NA_character_,
                 doseRx=doseRx, doseRxUnit=doseRxUnit,
                 isoDoseRx=isoDoseRx)
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
