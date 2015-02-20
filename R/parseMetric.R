## parse metric character strings into components and ensure validity
## optionally convert dose unit and volume unit
parseMetric <- function(x, doseUnit=NULL, volUnit=NULL) {
    ## remove whitespace and convert to upper case
    x <- toupper(removeWS(x))
    
    ## regular expression for components of a metric character string
    ## allow for V10%_CC or D10cc_% pattern for returning absolute volumes / relative doses
    pattern <- "^([DV])([.[:digit:]]+|MAX|MIN|MEAN|RX|SD)([%]|GY|CGY|CC)*(_[%]|_GY|_CGY|_CC)*.*"
    
    ## extract components
    DV      <- sub(pattern, "\\1", x)   # get a volume or a dose?
    valRef  <- sub(pattern, "\\2", x)   # given value at which to evaluate volume or dose
    unitRef <- sub(pattern, "\\3", x)   # measurement unit for given volume or dose
    unitDV_ <- sub(pattern, "\\4", x)   # measurement unit for output volume or dose

    ## remove _ from unitDV_
    unitDV <- ifelse(unitDV_ != "", sub("_", "", unitDV_), NA_character_)

    ## special dose cases: DMEAN, DSD, DMIN, DMAX, DRX
    specDose  <- valRef %in% c("MAX", "MIN", "MEAN", "RX", "SD")
    valRefNum <- ifelse((DV == "D") & specDose, NA_real_,      suppressWarnings(as.numeric(valRef)))
    unitRef   <- ifelse((DV == "D") & specDose, NA_character_, unitRef)

    ## convert absolute dose units if requested
    if(!is.null(doseUnit)) {
        doseUnit <- toupper(removeWS(doseUnit))

        ## output is dose -> reference is volume, just set unitDV
        ## don't convert relative volume
        idxD   <- (DV == "D") & (unitDV != "%")
        unitDV <- ifelse(idxD, doseUnit, unitDV)

        ## output is volume -> reference is dose
        ## don't convert relative dose
        idxV      <- (DV == "V") & (unitRef != "%")
        valRefNum <- ifelse(idxV,
                            valRefNum*suppressWarnings(getConvFac(paste0(unitRef, "2", doseUnit))),
                            valRefNum)

        valRef  <- ifelse(idxV, as.character(valRefNum), valRef)
        unitRef <- ifelse(idxV, doseUnit, unitRef)
    }
    
    ## convert absolute volume units if requested
    if(!is.null(volUnit)) {
        volUnit <- toupper(removeWS(volUnit))

        ## output is dose -> reference is volume
        ## don't convert relative volume
        ## consider special dose cases
        idxD      <- (DV == "D") & (unitRef != "%") & !specDose
        valRefNum <- ifelse(idxD,
                            valRefNum*suppressWarnings(getConvFac(paste0(unitRef, "2", volUnit))),
                            valRefNum)

        valRef  <- ifelse(idxD, as.character(valRefNum), valRef)
        unitRef <- ifelse(idxD, volUnit, unitRef)

        ## output is volume -> reference is dose, just set unitDV
        ## don't convert relative dose
        idxV   <- (DV == "V") & (unitDV != "%")
        unitDV <- ifelse(idxV, volUnit, unitDV)
    }

    ## canonical metric string
    unitDVStr  <- ifelse(is.na(unitDV),  "", paste0("_", unitDV))
    unitRefStr <- ifelse(is.na(unitRef), "", unitRef)
    metric     <- paste0(DV, valRef, unitRefStr, unitDVStr)

    ## check validity
    ## consider special cases DMEAN, DSD, DMIN, DMAX, DRX
    ## V -> %, Gy, cGy
    ## D -> %, CC
    validPattern <- grepl(pattern, x)
    validUnitRef <- ((DV == "D") &
                     ((unitRef %in% c("%", "CC")) |
                      (valRef  %in% c("MAX", "MIN", "MEAN", "RX", "SD")))) |
                    ((DV == "V") & (unitRef %in% c("%", "GY", "CGY")))
    validUnitDV  <- is.na(unitDV) |
                    ((DV == "D") & (unitDV  %in% c("%", "GY", "CGY"))) |
                    ((DV == "V") & (unitDV  %in% c("%", "CC")))

    valid <- validPattern & validUnitRef & validUnitDV
    if(!all(valid)) {
        warning(c("Pattern ", paste(x[!valid], collapse=", "), " is invalid"))
    }
 
    return(data.frame(metric=metric, valid=valid,
                      DV=DV, unitDV=unitDV,
                      valRef=valRef, valRefNum=valRefNum, unitRef=unitRef,
                      stringsAsFactors=FALSE))
}