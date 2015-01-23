## distance of point to linear interpolation of DVH curve
dvhDistance <-
function(x, DV) {
    ## relevant dose + volume data as matrix
    oneDistance <- function(x, oneDV) {
        if(any(is.na(oneDV))) {
            return(list(dstMin=NA_real_, ptMinD=NA_real_, ptMinV=NA_real_))
        }

        DVmat <- if(oneDV$volRel &&  oneDV$doseRel) {
            data.matrix(x$dvh[ , c("doseRel", "volumeRel")])
        } else if(  oneDV$volRel && !oneDV$doseRel) {
            data.matrix(x$dvh[ , c("dose",    "volumeRel")])
        } else if( !oneDV$volRel &&  oneDV$doseRel) {
            data.matrix(x$dvh[ , c("doseRel", "volume")])
        } else if( !oneDV$volRel && !oneDV$doseRel) {
            data.matrix(x$dvh[ , c("dose",    "volume")])
        }

        ## constraint point
        DVcpt <- c(oneDV$D, oneDV$V)

        #sp::spDistsN1(DVmat, DVcpt, longlat=FALSE)
        #DVlines   <- sp::Line(DVmat)
        #DVlinesL  <- sp::Lines(list(DVlines), ID="A")
        #DVlinesSL <- sp::SpatialLines(list(DVlinesL))
        #DVpt      <- sp::SpatialPoints(matrix(DVcpt, nrow=1))
        #rgeos::gDistance(DVlinesSL, DVpt)
        ## orthogonal projections of DV onto subspaces a'x
        ## defined by segments of x (DVcpt is treated as column vector)
        DVbase <- diff(DVmat)                 # basis vectors of segments
        lens   <- sqrt(rowSums(DVbase^2))     # length of basis vectors
        DVbu   <- diag(1/lens) %*% DVbase     # scaled to unit length
        op     <- DVbu %*% DVcpt              # orth proj in subspace coords

        ## check if projection is outside of segment -> not in [0, 1]
        outside <- (op < 0) | (op > 1)

        ## projections in DV coords aa'x
        DVbuL <- split(DVbu, f=seq_len(nrow(DVbu)))
        opOCL <- Map(function(v1, v2) { matrix(v1, ncol=1) %*% v2 }, DVbuL, op)

        ## add original coords back to move away from origin
        opOC   <- t(do.call("cbind", opOCL)) + DVmat[-nrow(DVmat), ]
        opOCin <- opOC[!outside, , drop=FALSE]

        ## distance of each point in DV to
        ## all points in DVmat as well as all non-outside orthogonal projections
        dstSupport <- sqrt(rowSums(sweep(DVmat,  2, DVcpt, "-")^2))
        dstProj    <- sqrt(rowSums(sweep(opOCin, 2, DVcpt, "-")^2))
        
        ## catch pathological cases
        if(!is.finite(min(dstSupport)) && !is.finite(min(dstProj))) {
            dstMin <- NA_real_
            ptMin  <- c(NA_real_, NA_real_)
        } else if(min(dstSupport) < suppressWarnings(min(dstProj))) {
            dstMin <- min(c(dstSupport, dstProj))
            ptMin  <- DVmat[which.min(dstSupport), ]
        } else {
            dstMin <- min(c(dstSupport, dstProj))
            ptMin  <- opOCin[which.min(dstProj), ]
        }

        return(list(dstMin=dstMin, ptMinD=ptMin[1], ptMinV=ptMin[2]))
    }

    ## get distance for each DV constraint point
    Map(oneDistance, list(x), split(DV, f=seq_len(nrow(DV))))
}
