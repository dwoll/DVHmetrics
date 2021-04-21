#####---------------------------------------------------------------------------
## parse character vector from Eclipse DVH file
parseEclipse <- function(x, planInfo=FALSE, courseAsID=FALSE, ...) {
    planInfo    <- as.character(planInfo)
    dots        <- list(...)
    uncertainty <- hasName(dots, "uncertainty") && dots$uncertainty

    ## function to extract one information element from a number of lines
    ## make sure only first : is matched -> not greedy
    getElem <- function(pattern, ll, trim=TRUE, iCase=FALSE, collWS=TRUE) {
        line <- ll[grep(pattern, ll)]
        elem <- sub("^.+?:[[:blank:]]+([[:alnum:][:punct:][:blank:]]+$)", "\\1",
                    line, ignore.case=iCase, perl=TRUE)
        elem <- if(trim) {
            trimws(elem, which="both")
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
        line <- trimws(ll[grep(pattern, ll)], which="both")
        pat  <- "^.+?:[^[:digit:]]*([[:digit:][:punct:]]+)(Gy|cGy|%)*([[:alnum:][:punct:][:blank:]]*)$"
        elem <- gsub(pat, "\\1", line, perl=TRUE, ignore.case=TRUE)
        num  <- if(length(elem) == 1L) {
            grp2 <- gsub(pat, "\\2", line, perl=TRUE, ignore.case=TRUE)
            grp3 <- gsub(pat, "\\3", line, perl=TRUE, ignore.case=TRUE)
            if(nzchar(grp2) || nzchar(grp3)) {
                header  <- x[seq_len(sStart[1]-1)]                        # header
                patID   <- getElem("^Patient ID[[:blank:]]*:",  header)   # patient id
                warning(paste(patID, ": Non-standard dose line found"))
            }
            
            as.numeric(trimws(elem))
        } else {
            NA_real_
        }
        
        if(is.na(num)) {
            header  <- x[seq_len(sStart[1]-1)]                        # header
            patID   <- getElem("^Patient ID[[:blank:]]*:",  header)   # patient id
            warning(paste(patID, ": No dose found"))
            return(num)
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
        line <- ll[grep("^(Prescribed|Total) dose.+:", ll)]
        elem <- sub("^.+\\[(GY|CGY)\\][[:blank:]]*:.+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        toupper(trimws(elem))
    }

    getVolUnit <- function(ll) {
        line <- ll[grep("^Volume.+:", ll)]
        if(length(line) >= 1L) {
            elem <- sub("^.+\\[(CM.+)\\][[:blank:]]*:.+", "\\1", line, perl=TRUE, ignore.case=TRUE)
            toupper(trimws(elem))
        } else {
            NA_character_
        }
    }

    getDVHtype <- function(ll) {
        line <- ll[grep("^Type.+:", ll)]
        elem <- sub("^Type.+:[[:blank:]]+(Cumulative|Differential).+", "\\1", line, perl=TRUE, ignore.case=TRUE)
        tolower(trimws(elem))
    }

    # "%A %d %B %Y %H:%M:%S"
    getDVHdate <- function(ll) {
        line <- ll[grep("^Date[[:blank:]]+:", ll)]
        elem <- trimws(sub("^Date[[:blank:]]+: (.+)$", "\\1", line))
        d <- if(grepl("^[[:digit:]]{2}\\.", elem)) {
            ## date format up to ARIA 12
            as.Date(substr(elem, start=1, stop=10), format="%d.%m.%Y")
        } else {
            ## date format ARIA 13+ - locale dependent
            ## remove ,
            elem <- gsub("[,.]", "", elem)
            as.Date(strptime(elem, format="%A %d %B %Y %H:%M:%S"))
        }
        
        if(is.na(d)) {
            elem
        } else {
            d
        }
    }

    if(uncertainty) {
        sec_indicator_string <- "(^Uncertainty plan: )|(^Structure: )"
        ## not quite start of sections because for uncertainty plans,
        ## both - "uncertainty plan" and "structure" appear in one block
        ## -> remove hits without preceding blank line
        sStart_tmp      <- grep(sec_indicator_string, x)
        sStart_tmp_prev <- sStart_tmp-1
        idx_uncertainty <- nzchar(x[sStart_tmp_prev])
        sStart          <- sStart_tmp[!idx_uncertainty]
        ## identify the uncertainty structures
        str_uncertainty <- which(idx_uncertainty)  - seq_len(sum(idx_uncertainty))
        str_regular     <- setdiff(seq_along(sStart), str_uncertainty)
    } else {
        sec_indicator_string <- "^Structure: [[:alnum:][:punct:]]"
        sStart <- grep(sec_indicator_string, x)     # start of sections
    }

    sLen <- diff(c(sStart, length(x)+1))      # length of sections
    if((length(sLen) < 1L) || all(sLen < 1L)) {
        stop("No structures found")
    }

    structList <- Map(function(start, len) x[start:(start+len-1)], sStart, sLen)

    ## extract file header and header info
    header  <- x[seq_len(sStart[1]-1)]                        # header
    patName <- getElem("Patient Name[[:blank:]]*:", header)   # patient name
    patID   <- getElem("^Patient ID[[:blank:]]*:",  header)   # patient id
    plan    <- getElem("^Plan[[:blank:]]?(sum)?:",  header)   # treatment plan
    course  <- getElem("^Course[[:blank:]]*:",      header)   # course

    ## generate new ID if courseAsID is set
    if(courseAsID) {
        patID <- removeWS(paste(patID, course, sep="_"))
    }

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

    DVHdate    <- getDVHdate(header)
    DVHtype    <- getDVHtype(header)
    isoDoseRx0 <- getElem("^% for dose \\(%\\):", header)
    ## check if sum plan
    isoDoseRx  <- if((length(isoDoseRx0) > 0) && (isoDoseRx0 != "not defined")) {
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

    doseUnit <- getDoseUnit(header)
    if(!grepl("^(GY|CGY)$", doseUnit)) {
        warning("Could not determine dose measurement unit")
        doseUnit <- NA_character_
    }

    doseRx0 <- getElem("^(Prescribed|Total) dose.*:", header)
    ## check if sum plan
    doseRx  <- if((length(doseRx0) > 0) && (doseRx0 != "not defined")) {
        doseRxUnit <- doseUnit
        getDose("^(Prescribed|Total) dose.*:", header)
    } else {                                        # sum plan
        ## doseRx is encoded in plan name
        if(tolower(planInfo) == "doserx") {
            doseRxUnit <- toupper(sub("^.+_([.[:digit:]]+)(GY|CGY)_.*", "\\2",
                                      plan, perl=TRUE, ignore.case=TRUE))
    
            if(!grepl("^(GY|CGY)$", doseRxUnit)) {
                warning("Could not determine dose Rx unit")
                doseRxUnit <- NA_character_
            }

            drx <- sub("^.+_([.[:digit:]]+)(GY|CGY)_.*", "\\1",
                       plan, perl=TRUE, ignore.case=TRUE)
            as.numeric(drx)
        } else {
            doseRxUnit <- NA_character_
            warning("No info on prescribed dose")
            NA_real_
        }
    }

    ## extract DVH from one structure section and store in a list
    ## with DVH itself as a matrix
    getDVH <- function(strct, info, is_uncertainty) {
        ## extract information from info list
        doseUnit   <- info$doseUnit
        doseRx     <- info$doseRx
        doseRxUnit <- info$doseRxUnit
        isoDoseRx  <- info$isoDoseRx

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
        if(is_uncertainty) {
            struct_org <- getElem("^Structure.*:", strct)
            structVol  <- info$str_volumes[[struct_org]]
            volumeUnit <- info$str_volUnits[[struct_org]]
            plan_uncertainty <- gsub("^Uncertainty plan: (.+) \\(variation .+\\)$", "\\1", strct[1])
            structure  <- paste(struct_org, trimws(plan_uncertainty), sep="_")
            
            doseMin   <- NA_real_
            doseMax   <- NA_real_
            doseAvg   <- NA_real_
            doseMed   <- NA_real_
            doseMode  <- NA_real_
            doseSD    <- NA_real_
        } else {
            structure <- getElem("^Structure.*:", strct)
            structVol <- as.numeric(getElem("^Volume.*:", strct))
            
            doseMin   <- getDose("^Min Dose.*:",    strct, doseRx)
            doseMax   <- getDose("^Max Dose.*:",    strct, doseRx)
            doseAvg   <- getDose("^Mean Dose.*:",   strct, doseRx)
            doseMed   <- getDose("^Median Dose.*:", strct, doseRx)
            doseMode  <- getDose("^Modal Dose.*:",  strct, doseRx)
            doseSD    <- getDose("^STD.*:",         strct, doseRx)
            
            volumeUnit <- getVolUnit(strct)
            volumeUnit <- if(grepl("^CM.+", volumeUnit, useBytes=TRUE)) {
                "CC"
            } else if(grepl("^%", volumeUnit)) {
                "PERCENT"
            } else {
                warning("Could not determine volume measurement unit")
                NA_character_
            }
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
                        split="\\[[[:alpha:]%]+\\]", fixed=FALSE, perl=TRUE))

        ## remove leading and trailing white space
        vars2 <- tolower(trimws(vars1))

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

        ## convert differential DVH (per unit dose) to cumulative
        ## and add differential DVH separately
        if(info$DVHtype == "differential") {
            DVH$dvh     <- convertDVH(dvh, toType="cumulative",
                                      toDoseUnit="asis", perDose=TRUE)
            DVH$dvhDiff <- dvh
        }

        ## set class
        class(DVH) <- "DVHs"
        return(DVH)
    }

    ## list of DVH data frames with component name = structure
    info <- list(patID=patID,
                 patName=patName,
                 date=DVHdate,
                 DVHtype=DVHtype,
                 plan=plan,
                 course=course,
                 quadrant=quadrant,
                 doseRx=doseRx,
                 doseRxUnit=doseRxUnit,
                 isoDoseRx=isoDoseRx,
                 doseUnit=doseUnit)

    ## if there are uncertainty structures, read corresponding normal
    ## structures first, extract volume, then read uncertainty structures
    dvhL <- if(uncertainty && (length(str_uncertainty) >= 1L)) {
        structList_regular     <- structList[str_regular]
        structList_uncertainty <- structList[str_uncertainty]
        
        dvhL_regular <- lapply(structList_regular,  getDVH, info=info,
                               is_uncertainty=FALSE)
        
        info$str_volumes <- setNames(Map(function(x) { x$structVol }, dvhL_regular),
                                     Map(function(x) { x$structure }, dvhL_regular))
        
        info$str_volUnits <- setNames(Map(function(x) { x$volumeUnit }, dvhL_regular),
                                      Map(function(x) { x$structure },  dvhL_regular))
        
        dvhL_uncertainty <- lapply(structList_uncertainty, getDVH, info=info,
                                   is_uncertainty=TRUE)
        
        c(dvhL_regular, dvhL_uncertainty)
        
    } else {
        lapply(structList, getDVH, info=info, is_uncertainty=FALSE)
    }
    dvhL <- Filter(Negate(is.null), dvhL)
    names(dvhL) <- sapply(dvhL, function(y) y$structure)
    if(length(unique(names(dvhL))) < length(dvhL)) {
        warning("Some structures have the same name - this can lead to problems")
    }

    class(dvhL) <- "DVHLst"
    attr(dvhL, which="byPat") <- TRUE

    return(dvhL)
}
