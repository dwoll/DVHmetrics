## TODO
## readDVH() -> if x is a directory, read all files in it
## differential -> cumulative: doses are interval mid-points
## make getMetrics(), showDVH(), checkConstraints(), showConstraints() work with RadOnc objects
## read Cadplan relVol, Pinnacle, Helax, TomoTherapy Hi-Art
## constraint wedges as custom geoms
## harmonize structures by reading equivalence file, possibly with regex

## build
#"c:\program files\r\r-3.1.0\bin\x64\Rcmd.exe" INSTALL --build DVHmetrics_0.1.tar.gz
#"c:\program files\r\r-3.1.0\bin\x64\Rcmd.exe" build DVHmetrics --resave-data --compact-vignettes="both"
#"c:\program files\r\r-3.1.0\bin\x64\Rcmd.exe" check DVHmetrics_0.1.tar.gz --as-cran
#install.packages("d:/daniel_work/workspace/DVHmetrics_0.1.tar.gz", repos=NULL, type="source")

#"c:\program files\r\r-devel\bin\x64\Rcmd.exe" INSTALL --build DVHmetrics_0.1.tar.gz
#"c:\program files\r\r-devel\bin\x64\Rcmd.exe" build DVHmetrics --resave-data --compact-vignettes="both"
#"c:\program files\r\r-devel\bin\x64\Rcmd.exe" check DVHmetrics_0.1.tar.gz --as-cran
#install.packages("h:/workspace/DVHmetrics_0.1.tar.gz", repos=NULL, type="source")

#####---------------------------------------------------------------------------
## trim whitespace on beginning/end of string
trimWS <- function(x, side="both")  {
    side <- match.arg(side, c("left", "right", "both"))
    pattern <- switch(side, left="^\\s+", right="\\s+$", both="^\\s+|\\s+$")
    gsub(pattern, "", x)
}

#####---------------------------------------------------------------------------
## collapse whitespace into one space
collWS <- function(x)  {
    gsub("[[:blank:]]+", " ", x)
}

#####---------------------------------------------------------------------------
## remove whitespace everywhere
removeWS <- function(x)  {
    gsub("[[:blank:]]+", "", x)
}

#####---------------------------------------------------------------------------
## parse character vector from Eclipse DVH file
parseEclipse <- function(x, planInfo=FALSE) {
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
        line <- ll[grep(pattern, ll)]
        elem <- sub("^.+?:[[:blank:]]+([[:alnum:][:punct:]]+[[:blank:]]*$)", "\\1", line, perl=TRUE)
        num  <- trimWS(elem)
        if(grepl("\\[%\\]", line)) {
            ## relative dose
            doseRx * as.numeric(num)/100
        } else {
            as.numeric(num)
        }
    }

    getDoseUnit <- function(ll) {
        line <- ll[grep("^Prescribed dose.+:", ll)]
        elem <- sub("^.+\\[(GY|CGY)\\][[:blank:]]*:.+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        toupper(trimWS(elem))
    }

    getVolUnit <- function(ll) {
        line <- ll[grep("^Volume.+:", ll)]
        elem <- sub("^.+\\[(CM.+)\\][[:blank:]]*:.+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        toupper(trimWS(elem))
    }

    getDVHtype <- function(ll) {
        line <- ll[grep("^Type.+:", ll)]
        elem <- sub("^Type.+:[[:blank:]]+(Cumulative|Differential).+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        tolower(trimWS(elem))
    }

    sStart <- grep("^Structure: [[:alnum:]]", x)   # start of sections
    sLen   <- diff(c(sStart, length(x)+1))         # length of sections
    if((length(sLen) < 1L) || all(sLen < 1L)) {
        stop("No structures found")
    }

    structList <- Map(function(start, len) x[start:(start+len-1)], sStart, sLen)

    ## extract file header and header info
    header  <- x[seq_len(sStart[1]-1)]                        # header
    patName <- getElem("Patient Name[[:blank:]]*:", header)   # patient name
    patID   <- getElem("^Patient ID[[:blank:]]*:",  header)   # patient id
    plan    <- getElem("^Plan[[:blank:]]?(sum)?:",  header)   # treatment plan

    ## some information specific to PASSOS - quadrant is coded in plan
    quadrant <- if(planInfo == "QuadrantUlm") {
        ## quadrant from plan
        sub("^Ma[[:alpha:]]*(Li|Re)[[:blank:]]*([[:alpha:]]+)", "\\2", plan, perl=TRUE)
    } else if(planInfo == "QuadrantMainz") {
        patThxw <- "^PA_1-2_(L|R)Thoraxwa"   # Thoraxwand pattern
        if(grepl(patThxw, plan, ignore.case=TRUE, perl=TRUE)) {
            "T"
        } else {
            ## possible quadrants
            quadsHave <- c(1:4, "O", "U", "L", "M", "Z")
            quads     <- c(quadsHave, "-", "X") # with unknown
            patQuad   <- "^PA_1-2_R(-|[[:alnum:]]+)_L(-|[[:alnum:]]+)"
            quadR     <- sub(patQuad, "\\1", plan, ignore.case=TRUE, perl=TRUE)
            quadL     <- sub(patQuad, "\\2", plan, ignore.case=TRUE, perl=TRUE)
            if(!all(toupper(quadR) %in% quads)) {
                quadR <- "X"
            }

            if(!all(toupper(quadL) %in% quads)) {
                quadL <- "X"
            }

            if(quadR == "-")  {
                quadL
            } else if(quadL == "-") {
                quadR
            } else if((quadR %in% quadsHave) && (quadL %in% quadsHave)) {
                "B"    # both sides
            } else {
                "X"    # unknown
            }
        }
    } else {
        NA_character_
    }

    DVHtype    <- getDVHtype(header)
    isoDoseRx0 <- getElem("^% for dose \\(%\\):", header)
    isoDoseRx  <- if(isoDoseRx0 != "not defined") { # check if sum plan
        as.numeric(isoDoseRx0)
    } else {                                        # sum plan -> use plan info?
        if(tolower(planInfo) == "doserx") {
            warning("Iso-dose-Rx is assumed to be 100")
            100
        } else {
            warning("No info on % for dose")
            NA_real_
        }
    }

    doseRx0 <- getElem("^Prescribed dose.*:", header)
    doseRx  <- if(doseRx0 != "not defined") {       # check if sum plan
        getDose("^Prescribed dose.*:", header)
    } else {                                        # sum plan
        ## doseRx is encoded in plan name
        if(tolower(planInfo) == "doserx") {
            drx <- sub("^[[:alnum:]]+_([.[:digit:]]+)(GY|CGY)_[[:alnum:]]*", "\\1",
                       plan, perl=TRUE, ignore.case=TRUE)
            as.numeric(drx)
        } else {
            warning("No info on prescribed dose")
            NA_real_
        }
    }

    doseUnit <- getDoseUnit(header)
    if(!grepl("^(GY|CGY)$", doseUnit)) {
        warning("Could not determine dose measurement unit")
        doseUnit <- NA_character_
    }

    ## extract DVH from one structure section and store in a list
    ## with DVH itself as a data frame
    getDVH <- function(strct, info) {
        ## extract information from info list
        plan      <- info$plan
        quadrant  <- info$quadrant
        doseRx    <- info$doseRx
        isoDoseRx <- info$isoDoseRx

        ## extract structure, volume, dose min, max, mean, median and sd
        structure <- getElem("^Structure.*:", strct)
        structVol <- as.numeric(getElem("^Volume.*:", strct))
        doseMin   <- getDose("^Min Dose.*:",    strct, doseRx)
        doseMax   <- getDose("^Max Dose.*:",    strct, doseRx)
        doseAvg   <- getDose("^Mean Dose.*:",   strct, doseRx)
        doseMed   <- getDose("^Median Dose.*:", strct, doseRx)
        doseMode  <- getDose("^Modal Dose.*:",  strct, doseRx)
        doseSD    <- getDose("^STD.*:",         strct, doseRx)

        volumeUnit <- getVolUnit(strct)
        volumeUnit <- if(grepl("^CM.+", volumeUnit)) {
            "CC"
        } else {
            warning("Could not determine volume measurement unit")
            NA_character_
        }

        ## find DVH
        ## DVH column headers
        colHead  <- grep("DOSE \\[(%|GY|CGY)\\].+VOLUME", strct,
                         ignore.case=TRUE, perl=TRUE)
        dvhStart <- colHead+1                 # first numeric line of DVH
        dvhLen   <- length(strct) - dvhStart
        if((length(dvhLen) < 1L) || dvhLen < 1L) {
            stop("No DVH data found")
        }

        ## column headers
        vars1 <- unlist(strsplit(strct[colHead],
                        split="\\[[[:alpha:]%]+\\]", fixed=FALSE, perl=TRUE))

        ## remove leading and trailing white space
        vars2 <- tolower(trimWS(vars1))

        ## make sure we recognize all columns in the DVH
        patDose    <- "^dose"
        patDoseRel <- "^relative dose"
        patVol     <- "^structure volume"
        patVolRel  <- "^ratio[[:alpha:][:blank:]]+volume"
        patVolD    <- "^dvolume / ddose"  # differential DVH
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
        vars3[grep(patVolD,    vars2)] <- "volume"

        ## extract DVH as a matrix and store preceding information
        ## read line length(strct) for cases where file does not end with a
        ## blank line -> this will then be last DVH line, otherwise blank
        con <- textConnection(strct[dvhStart:length(strct)])
        dvh <- data.matrix(read.table(con,
                                      header=FALSE, stringsAsFactors=FALSE,
                                      colClasses=rep("numeric", length(vars3)),
                                      comment.char="", nrows=dvhLen))
        close(con)

        ## rename
        colnames(dvh) <- vars3

        ## catch special case: structVol is 0.0 due to limited precision
        if("volume" %in% vars3) {
            structVol <- if(info$DVHtype == "cumulative") {
                max(c(structVol, dvh[ , "volume"]))
            } else {
                ## reconstruct volumes -> volume is per gray -> mult with bin width
                volBin <- dvh[ , "volume"]*diff(c(-dvh[1, "dose"], dvh[ , "dose"]))
                max(c(structVol, sum(volBin)))
            }
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
            dvh <- cbind(dvh, dose=dvh[ , "doseRel"]*doseRx / isoDoseRx)
            # (doseRx/(isoDoseRx/100))*(dvh$doseRel/100)
        } else if(("dose" %in% vars3) && !("doseRel" %in% vars3)) {
            dvh <- cbind(dvh, doseRel=dvh[ , "dose"]*isoDoseRx / doseRx)
            # 100*(dvh$dose/(doseRx/(isoDoseRx/100)))
        }

        DVH <- list(dvh=dvh,
                    patName=info$patName,
                    patID=info$patID,
                    DVHtype=info$DVHtype,
                    plan=plan,
                    quadrant=quadrant,
                    structure=structure,
                    structVol=structVol,
                    doseUnit=info$doseUnit,
                    volumeUnit=volumeUnit,
                    doseMin=doseMin,
                    doseMax=doseMax,
                    doseRx=doseRx,
                    isoDoseRx=isoDoseRx,
                    doseAvg=doseAvg,
                    doseMed=doseMed,
                    doseMode=doseMode,
                    doseSD=doseSD)

        ## convert differential DVH to cumulative
        ## and add differential DVH separately
        if(info$DVHtype == "differential") {
            DVH$dvh     <- dvhConvert(dvh, toType="cumulative", toDoseUnit="asis")
            DVH$dvhDiff <- dvh
        }

        ## set class
        class(DVH) <- c("DVHs", "list")
        return(DVH)
    }

    ## list of DVH data frames with component name = structure
    info <- list(patID=patID, patName=patName,
                 DVHtype=DVHtype, plan=plan, quadrant=quadrant,
                 doseRx=doseRx, isoDoseRx=isoDoseRx, doseUnit=doseUnit)
    dvhL <- lapply(structList, getDVH, info=info)
    names(dvhL) <- sapply(dvhL, function(y) y$structure)
    class(dvhL) <- c("DVHLst", "list")
    attr(dvhL, which="byPat") <- TRUE

    return(dvhL)
}

## parse character vector from Cadplan DVH file
parseCadplan <- function(x, planInfo=FALSE) {
    ## function to extract one information element from a number of lines
    getElem <- function(pattern, ll, trim=TRUE, iCase=FALSE) {
        line <- ll[grep(pattern, ll)]
        elem <- sub("^.+?:[[:blank:]]*([[:alnum:][:punct:][:blank:]]+$)", "\\1",
                    line, ignore.case=iCase, perl=TRUE)
        if(trim) {
            trimWS(elem, side="both")
        } else {
            elem
        }
    }

    getPlan <- function(pattern, ll, trim=TRUE, iCase=FALSE, collWS=TRUE) {
        line <- ll[grep(pattern, ll)]
        elem <- sub("^[[:alpha:]]+[[:blank:]]+([[:alnum:]]+[[:alnum:][:blank:]]*)", "\\1",
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

    getDose <- function(pattern, ll, doseRx, percent=TRUE) {
        line <- ll[grep(pattern, ll)]
        elem <- sub("^.+?:[[:blank:]]+([[:alnum:][:punct:]]+[[:blank:]]*$)", "\\1", line, perl=TRUE)
        num  <- trimWS(elem)
        if(percent && any(grepl("%", line))) {
            doseRx * as.numeric(num)/100
        } else {
            as.numeric(num)
        }
    }

    getDoseUnit <- function(ll) {
        line <- ll[grep("^Prescr\\. dose.+:", ll)]
        elem <- sub("^.+\\((GY|CGY)\\)[[:blank:]]*:.+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        toupper(trimWS(elem))
    }

    getVolUnit <- function(ll) {
        line <- ll[grep("^Volume.+:", ll)]
        elem <- sub("^.+\\((CC)\\)[[:blank:]]*:.+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        toupper(trimWS(elem))
    }

    getDVHtype <- function(ll) {
        line <- ll[grep("^(Cumulative|Differential) Dose Volume.+", ll)]
        elem <- sub("^(Cumulative|Differential).+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        tolower(trimWS(elem))
    }

    ## split file into list of structure sections
    sStart <- grep("^Histogram[[:blank:]]*:[[:blank:]]*[[:alnum:]]+", x)  # start of sections
    sLen   <- diff(c(sStart, length(x)+1))        # length of sections
    if((length(sLen) < 1L) || all(sLen < 1L)) {
        stop("No structures found")
    }

    structList <- Map(function(start, len) x[start:(start+len-1)], sStart, sLen)

    ## extract file header and header info
    header  <- x[seq_len(sStart[1]-1)]                       # header
    patName <- getElem("Patient Name[[:blank:]]*:", header)  # patient name
    patID   <- getElem("^Patient ID[[:blank:]]*:",  header)  # patient id
    plan    <- getPlan("^PLAN", header, iCase=TRUE, collWS=TRUE)  # treatment plan
    DVHtype <- getDVHtype(header)

    ## extract DVH from one structure section and store in a list
    ## with DVH itself as a data frame
    getDVH <- function(strct, info) {
        plan <- info$plan

        ## extract structure, prescribed dose, volume, completed to 100% dose Rx,
        ## dose min, max, mean, median, mode, and sd
        structure <- getElem("^Histogram.*:", strct)

        isoDoseRx0 <- getElem("^% for dose[[:blank:]]*:", strct)
        isoDoseRx  <- if(isoDoseRx0 != "not defined") { # check if sum plan
            as.numeric(isoDoseRx0)
        } else {                                        # sum plan -> use plan info?
            if(tolower(planInfo) == "doserx") {
                warning("Iso-dose-Rx is assumed to be 100")
                100
            } else {
                warning("No info on % for dose")
                NA_real_
            }
        }

        doseRx0 <- getElem("^Prescr\\. dose.*:", strct)
        doseRx  <- if(doseRx0 != "not defined") {          # check if sum plan
            getDose("^Prescr\\. dose.*:", strct)
        } else {                                        # sum plan
            ## doseRx is encoded in plan name
            if(tolower(planInfo) == "doserx") {
                drx <- sub("^[[:alnum:]]+_([.[:digit:]]+)(GY|CGY)_[[:alnum:]]*", "\\1",
                           plan, perl=TRUE, ignore.case=TRUE)
                as.numeric(drx)
            } else {
                warning("No info on prescribed dose")
                NA_real_
            }
        }

        doseUnit <- getDoseUnit(strct)
        if(!grepl("^(GY|CGY)$", doseUnit)) {
            doseUnit <- NA_character_
            warning("Could not determine dose measurement unit")
        }

        volumeUnit <- getVolUnit(strct)
        if(!grepl("^CC$", volumeUnit)) {
            volumeUnit <- NA_character_
            warning("Could not determine volume measurement unit")
        }

        structVol <- as.numeric(getElem("^Volume.*:", strct))
        doseMin   <- getDose("^Dose minimum.*:", strct, doseRx)
        doseMax   <- getDose("^Dose maximum.*:", strct, doseRx)
        doseAvg   <- getDose("^Dose mean.*:",    strct, doseRx)
        doseMed   <- getDose("^Dose median.*:",  strct, doseRx)
        doseMod   <- getDose("^Dose modal.*:",   strct, doseRx)
        doseSD    <- getDose("^Standard dev.*:", strct, doseRx, percent=FALSE)

        ## find DVH
        ## DVH column headers
        colHead  <- grep("DOSE[[:blank:]]*\\((%|GY|CGY)\\).+VOLUME", strct,
                         ignore.case=TRUE, perl=TRUE)
        dvhStart <- colHead+1            # first numeric line of DVH
        dvhLen   <- length(strct) - dvhStart
        if((length(dvhLen) < 1L) || dvhLen < 1L) {
            stop("No DVH data found")
        }

        ## column headers
        vars1 <- unlist(strsplit(strct[colHead],
                        split="\\([[:alpha:]%]+\\)", fixed=FALSE, perl=TRUE))
        ## remove leading and trailing white space
        vars2 <- tolower(trimWS(vars1))

        ## make sure we recognize all columns in the DVH
        patDose    <- "^dose"
        patDoseRel <- "^relative dose"
        patVol     <- "^volume.+cm3"
        hits <- sum(c(grepl(patDose, vars2), grepl(patDoseRel, vars2), grepl(patVol, vars2)))
        if(length(vars2) != hits) {
        	stop(c("Could not identify all DVH columns"),
        		 paste(vars2, collapse=", "))
        }

        ## replace column headers
        vars3 <- vars2
        vars3[grep(patDose,    vars2)] <- "dose"
        vars3[grep(patDoseRel, vars2)] <- "doseRel"
        vars3[grep(patVol,     vars2)] <- "volume"

        ## extract DVH as a matrix and store preceding information
        ## read line length(strct) for cases where file does not end with a
        ## blank line -> this will then be last DVH line, otherwise blank
        con <- textConnection(strct[dvhStart:length(strct)])
        dvh <- data.matrix(read.table(con,
                                      header=FALSE, stringsAsFactors=FALSE,
                                      colClasses=rep("numeric", length(vars3)),
                                      comment.char="", nrows=dvhLen))
        close(con)

        ## rename
        colnames(dvh) <- vars3

        ## catch special case: structVol is 0.0 due to limited precision
        if("volume" %in% vars3) {
            structVol <- if(info$DVHtype == "cumulative") {
                max(c(structVol, dvh[ , "volume"]))
            } else {
                ## reconstruct volumes -> volume is per gray -> mult with bin width
                volBin <- dvh[ , "volume"]*diff(c(-dvh[1, "dose"], dvh[ , "dose"]))
                max(c(structVol, sum(volBin)))
            }
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
            dvh <- cbind(dvh, dose=dvh[ , "doseRel"]*doseRx / isoDoseRx)
            # (doseRx/(isoDoseRx/100))*(dvh$doseRel/100)
        } else if(("dose" %in% vars3) && !("doseRel" %in% vars3)) {
            dvh <- cbind(dvh, doseRel=dvh[ , "dose"]*isoDoseRx / doseRx)
            # 100*(dvh$dose/(doseRx/(isoDoseRx/100)))
        }

        DVH <- list(dvh=dvh,
                    patID=info$patID,
                    patName=info$patName,
                    DVHtype=info$DVHtype,
                    plan=plan,
                    structure=structure,
                    structVol=structVol,
                    doseUnit=doseUnit,
                    volumeUnit=volumeUnit,
                    doseMin=doseMin,
                    doseMax=doseMax,
                    doseRx=doseRx,
                    isoDoseRx=isoDoseRx,
                    doseAvg=doseAvg,
                    doseMed=doseMed,
                    doseSD=doseSD)

        ## convert differential DVH to cumulative
        ## and add differential DVH separately
        if(info$DVHtype == "differential") {
            DVH$dvh     <- dvhConvert(dvh, toType="cumulative", toDoseUnit="asis")
            DVH$dvhDiff <- dvh
        }

        ## set class
        class(DVH) <- c("DVHs", "list")
        return(DVH)
    }

    ## list of DVH data frames with component name = structure
    info <- list(patID=patID, patName=patName, plan=plan, DVHtype=DVHtype)
    dvhL <- lapply(structList, getDVH, info)
    names(dvhL) <- sapply(dvhL, function(y) y$structure)
    class(dvhL) <- c("DVHLst", "list")
    attr(dvhL, which="byPat") <- TRUE

    return(dvhL)
}

## returns a list (1 component per DVH file) of lists (1 component = 1 list per structure)
readDVH <- function(x, type=c("Eclipse", "Cadplan"), planInfo=FALSE) {
    type <- match.arg(type)

    dvhRawL <- if(missing(x)) {
        parseDVH()
    } else {
        parseDVH(x)
    }

    dvhLL <- if(type == "Eclipse") {
        lapply(dvhRawL, parseEclipse, planInfo=planInfo)
    } else if(type == "Cadplan") {
        lapply(dvhRawL, parseCadplan, planInfo=planInfo)
    }

    ## organized by patient (top level)
    attr(dvhLL, which="byPat") <- TRUE
    class(dvhLL) <- c("DVHLstLst", "list")

    return(dvhLL)
}
