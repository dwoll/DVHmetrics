get_name_pair <- function(x, y, sep=" <-> ") {
    paste0(x, sep, y) 
}

get_name_elem <- function(x, which=1L, sep=" <-> ") {
    stopifnot(which %in% c(1L, 2L))
    stopifnot(grepl(sep, x))
    if(which == 1) {
        gsub(paste0("^(.+)", sep, ".+$"), "\\1", x)
    } else {
        gsub(paste0("^.+", sep, "(.+$)"), "\\1", x)
    }
}

reconstruct_mesh <- function(x,
                             method=c("No", "AFS", "SSS", "Poisson"),
                             jetSmoothing=NULL,
                             scaleIterations=1,
                             neighbors=12,
                             samples=300,
                             separateShells=FALSE,
                             forceManifold=TRUE,
                             borderAngle=45,
                             normals="jet(12)",
                             spacing="ave(12)",
                             sm_angle=20,
                             sm_radius=30,
                             sm_distance=0.375) {
    method <- match.arg(tolower(method),
                        choices=c("no", "afs", "sss", "poisson"))
    
    ## check functions available in cgalMeshes version
    cgalMeshes_has     <- ls(getNamespace("cgalMeshes"), all.names=TRUE)
    cgalMeshes_version <- packageVersion("cgalMeshes")
    cgalMeshes_atleast110 <- compareVersion("1.1.0", as.character(cgalMeshes_version))
    
    if(method == "afs") {
        if(cgalMeshes_atleast110 <= 0) {
            AFSreconstruction(x$getVertices(),
                              jetSmoothing=jetSmoothing) # for cgalMeshes 1.1.0
        } else {
            AFSreconstruction(x$vertices()) # for cgalMeshes 1.0.0
        }
    } else if(method == "sss") {
        if("SSSreconstruction" %in% cgalMeshes_has) {
            cgalMeshes::SSSreconstruction(x$getVertices(),
                                          scaleIterations=scaleIterations,
                                          neighbors=neighbors,
                                          samples=samples,
                                          separateShells=separateShells,
                                          forceManifold=forceManifold,
                                          borderAngle=borderAngle)
        } else {
            warning("SSS reconstruction not implemented. Using AFS instead.")
            if(cgalMeshes_atleast110 <= 0) {
                AFSreconstruction(x$getVertices(),
                                  jetSmoothing=jetSmoothing) # for cgalMeshes 1.1.0
            } else {
                AFSreconstruction(x$vertices()) # for cgalMeshes 1.0.0
            }
        }
    } else if(method == "poisson") {
        if("PoissonReconstruction" %in% cgalMeshes_has) {
            cgalMeshes::PoissonReconstruction(x$getVertices(),
                                              normals=normals,
                                              spacing=spacing,
                                              sm_angle=sm_angle,
                                              sm_radius=sm_radius,
                                              sm_distance=sm_distance)
        } else {
            warning("Poisson reconstruction not implemented. Using AFS instead.")
            if(cgalMeshes_atleast110 <= 0) {
                AFSreconstruction(x$getVertices(),
                                  jetSmoothing=jetSmoothing) # for cgalMeshes 1.1.0
            } else {
                AFSreconstruction(x$vertices()) # for cgalMeshes 1.0.0
            }
        }
    } else {
        stop("Wrong reconstruction method")
    }
}

read_mesh_one <- function(x,
                          name,
                          fix_issues=TRUE,
                          reconstruct=c("No", "AFS", "SSS", "Poisson"),
                          jetSmoothing=NULL,
                          scaleIterations=1,
                          neighbors=12,
                          samples=300,
                          separateShells=FALSE,
                          forceManifold=TRUE,
                          borderAngle=45,
                          normals="jet(12)",
                          spacing="ave(12)",
                          sm_angle=20,
                          sm_radius=30,
                          sm_distance=0.375) {
    reconstruct <- match.arg(tolower(reconstruct),
                             choices=c("no", "afs", "sss", "poisson"))
    
    mesh_name <- if(missing(name)) {
        basename(file_path_sans_ext(x))
    } else {
        basename(file_path_sans_ext(name))
    }
    
    mesh <- cgalMesh$new(x)
    
    if(reconstruct != "no") {
        mesh_r <- reconstruct_mesh(mesh,
                                   method=reconstruct,
                                   jetSmoothing=jetSmoothing,
                                   neighbors=neighbors,
                                   samples=samples,
                                   separateShells=separateShells,
                                   forceManifold=forceManifold,
                                   borderAngle=borderAngle,
                                   normals=normals,
                                   spacing=spacing,
                                   sm_angle=sm_angle,
                                   sm_radius=sm_radius,
                                   sm_distance=sm_distance)
        mesh <- mesh_r
    }
    
    diag_nsi    <- !mesh$selfIntersects()
    diag_closed <- mesh$isClosed()
    diag_bv     <- mesh$boundsVolume()
    
    issues <- c("self intersects", "not closed", "does not bound volume")
    
    if(!all(diag_nsi, diag_closed, diag_bv)) {
        warn_str <- paste0("Mesh ", mesh_name, " has these issues: ",
                           paste(issues[!c(diag_nsi, diag_closed, diag_bv)],
                                 collapse=", "))
        
        if(fix_issues) {
            warn_str <- paste0(warn_str, ". Trying to fix.")
            warning(warn_str)
            
            if(!diag_closed) {
                if(reconstruct == "no") {
                    reconstruct <- "afs"
                }
                
                mesh_r <- reconstruct_mesh(mesh, method=reconstruct, spacing=spacing)
                mesh   <- mesh_r
            }
            
            if(!mesh$boundsVolume()) {
                mesh$orientToBoundVolume()
            }
            
            if(mesh$selfIntersects()) {
                mesh$removeSelfIntersections()
            }
        }
    }
       
    vol_0 <- try(mesh$volume())
    ctr_0 <- try(mesh$centroid())
    vol <- if(!inherits(vol_0, "try-error")) {
        abs(vol_0)
    } else {
        NA_real_
    }
    
    ctr <- if(!inherits(ctr_0, "try-error")) {
        ctr_0
    } else {
        rep(NA_real_, 3)
    }
    
    list(name    =mesh_name,
         mesh    =mesh,
         volume  =vol,
         centroid=ctr)
}

read_mesh_obs <- function(x,
                          name,
                          fix_issues=TRUE,
                          reconstruct=c("No", "AFS", "SSS", "Poisson"),
                          jetSmoothing=NULL,
                          scaleIterations=1,
                          neighbors=12,
                          samples=300,
                          separateShells=FALSE,
                          forceManifold=TRUE,
                          borderAngle=45,
                          normals="jet(12)",
                          spacing="ave(12)",
                          sm_angle=20,
                          sm_radius=30,
                          sm_distance=0.375) {
    reconstruct <- match.arg(tolower(reconstruct),
                             choices=c("no", "afs", "sss", "poisson"))
    
    mesh_names <- if(missing(name)) {
        basename(tools::file_path_sans_ext(x))
    } else {
        basename(tools::file_path_sans_ext(name))
    }
    
    meshL <- lapply(seq_along(x), function(i) {
        read_mesh_one(x[i],
                      name=mesh_names[i],
                      fix_issues=fix_issues,
                      reconstruct=reconstruct,
                      jetSmoothing=jetSmoothing,
                      neighbors=neighbors,
                      samples=samples,
                      separateShells=separateShells,
                      forceManifold=forceManifold,
                      borderAngle=borderAngle,
                      normals=normals,
                      spacing=spacing,
                      sm_angle=sm_angle,
                      sm_radius=sm_radius,
                      sm_distance=sm_distance)
    })
    
    setNames(meshL, mesh_names)
}

read_mesh <- function(x,
                      name,
                      fix_issues=TRUE,
                      reconstruct=c("No", "AFS", "SSS", "Poisson"),
                      jetSmoothing=NULL,
                      scaleIterations=1,
                      neighbors=12,
                      samples=300,
                      separateShells=FALSE,
                      forceManifold=TRUE,
                      borderAngle=45,
                      normals="jet(12)",
                      spacing="ave(12)",
                      sm_angle=20,
                      sm_radius=30,
                      sm_distance=0.375) {
    reconstruct <- match.arg(tolower(reconstruct),
                             choices=c("no", "afs", "sss", "poisson"))
    
    obs_names <- if(missing(name)) {
        names(x)
    } else {
        names(name)
    }
    
    mesh_names <- if(missing(name)) {
        lapply(x, function(y) { basename(tools::file_path_sans_ext(y)) })
    } else {
        lapply(x, function(y) { basename(tools::file_path_sans_ext(name)) })
    }
    
    if(is.null(jetSmoothing)) {
        jetSmoothing <- list(NULL)
    }
    
    Map(read_mesh_obs,
        setNames(x, obs_names),
        mesh_names,
        fix_issues=fix_issues,
        reconstruct=reconstruct,
        jetSmoothing=jetSmoothing,
        neighbors=neighbors,
        samples=samples,
        separateShells=separateShells,
        forceManifold=forceManifold,
        borderAngle=borderAngle,
        normals=normals,
        spacing=spacing,
        sm_angle=sm_angle,
        sm_radius=sm_radius,
        sm_distance=sm_distance)
    }

get_mesh_info_one <- function(x) {
    mesh_list <- x[["mesh"]]$getMesh(rgl=FALSE)
    data.frame(name=x[["name"]],
               n_verts=nrow(mesh_list[["vertices"]]),
               n_faces=ncol(mesh_list[["faces"]]),
               volume=x[["volume"]],
               ctr_x=x[["centroid"]][1],
               ctr_y=x[["centroid"]][2],
               ctr_z=x[["centroid"]][3])
}

get_mesh_info_obs <- function(x) {
    d_out <- do.call("rbind", Map(get_mesh_info_one, x))
    rownames(d_out) <- NULL
    d_out
}

get_mesh_info <- function(x) {
    x_lens <- lengths(x)
    d_out  <- do.call("rbind", Map(get_mesh_info_obs, x))
    rownames(d_out) <- NULL
    cbind(observer=rep(names(x), times=x_lens),
          d_out)
}

## starting from list of observers, each with a list of meshes
## generate all observer-pairs for each corresponding mesh-list entry
get_mesh_pairs <- function(x, sep=" <-> ", names_only=FALSE) {
    if(length(x) <= 1L) { stop("Need more than 1 mesh for comparisons") }
    
    ## number of meshes per observer
    n_obs_meshes <- lengths(x) # may be different per observer
    n_meshes     <- max(n_obs_meshes)
    
    ## for given pair, put corresponding meshes in a list
    get_pair_mesh <- function(pair, mesh) {
        idx1  <- pairs_idx[pair, 1] # index observer 1
        idx2  <- pairs_idx[pair, 2] # index observer 2
        obs1  <- x[[idx1]]          # observer 1
        obs2  <- x[[idx2]]          # observer 2
        ## do both observers have the mesh?
        if((length(obs1) >= mesh) && (length(obs2) >= mesh)) {
            mesh_1 <- x[[idx1]][[mesh]]
            mesh_2 <- x[[idx2]][[mesh]]
            
            ## add group information -> same structure
            if(names_only) {
                list(name=get_name_pair(mesh_1[["name"]], mesh_2[["name"]], sep=sep),
                     group=sprintf("strct_%.3d", mesh))
            } else {
                list(name=get_name_pair(mesh_1[["name"]], mesh_2[["name"]], sep=sep),
                     mesh_1=mesh_1,
                     mesh_2=mesh_2,
                     group=sprintf("strct_%.3d", mesh))
            }
        } else {
            NULL
        }
    }
    
    pairs_idx <- if(length(x) >= 2L) {
        t(combn(seq_along(x), 2))
    } else {
        matrix(c(1, 1), ncol=2)
    }
    
    ll_outer <- lapply(seq_len(n_meshes), function(idx_mesh) {
        lapply(seq_len(nrow(pairs_idx)), function(idx_pair) { get_pair_mesh(idx_pair, idx_mesh) })
    })
    
    ## weed out NULL components
    ll <- Filter(Negate(is.null), unlist(ll_outer, recursive=FALSE))
    pair_names <- lapply(ll, function(x) { x[["name"]] })
    setNames(ll, pair_names)
}

## union and intersection for list of two meshes x
get_mesh_ui_pair <- function(x, boov=FALSE) {
    if(!boov) {
        m_union     <- try(x[["mesh_1"]][["mesh"]]$union(       x[["mesh_2"]][["mesh"]]))
        m_intersect <- try(x[["mesh_1"]][["mesh"]]$intersection(x[["mesh_2"]][["mesh"]]))
        ui_ok       <- !(inherits(m_union,     "try-error") ||
                         inherits(m_intersect, "try-error"))
    } else {
        have_Boov <- requireNamespace("Boov", quietly=TRUE, partial=TRUE)
        if(!have_Boov) { warning("Package 'Boov' required for 'boov=TRUE' but not found.") }
        stopifnot(have_Boov)

        m1_rgl <- x[["mesh_1"]][["mesh"]]$getMesh(rgl=TRUE)
        m2_rgl <- x[["mesh_2"]][["mesh"]]$getMesh(rgl=TRUE)
        
        m_union_0     <- try(Boov::MeshesUnion(       list(m1_rgl, m2_rgl), clean=TRUE))
        m_intersect_0 <- try(Boov::MeshesIntersection(list(m1_rgl, m2_rgl), clean=TRUE))
        ui_ok         <- !(inherits(m_union_0,     "try-error") ||
                           inherits(m_intersect_0, "try-error"))
        
        if(ui_ok) {
            m_union_rgl     <- Boov::toRGL(m_union_0)
            m_intersect_rgl <- Boov::toRGL(m_intersect_0)
            
            m_union     <- try(cgalMesh$new(m_union_rgl))
            m_intersect <- try(cgalMesh$new(m_intersect_rgl))
            
            ## intersection might be empty
            ui_ok <- !(inherits(m_union,     "try-error") ||
                       inherits(m_intersect, "try-error"))
        }
    }
    
    if(!ui_ok) {
        m_union     <- NULL
        m_intersect <- NULL
        vol_u       <- NA_real_
        vol_i       <- NA_real_
    } else {
        if(m_union$selfIntersects()) {
            m_union$removeSelfIntersections()
        }
        
        if(m_intersect$selfIntersects()) {
            m_intersect$removeSelfIntersections()
        }
        
        if(!m_union$boundsVolume()) {
            m_union$orientToBoundVolume()
        }
        
        if(!m_intersect$boundsVolume()) {
            m_intersect$orientToBoundVolume()
        }
        
        vol_u_0 <- try(m_union$volume())
        vol_i_0 <- try(m_intersect$volume())
        
        vol_u <- if(inherits(vol_u_0, "try-error") || is.na(vol_u_0)) {
            NA_real_
        } else {
            vol_u_0
        }

        vol_i <- if(inherits(vol_i_0, "try-error") || is.na(vol_i_0)) {
            NA_real_
        } else {
            vol_i_0
        }
    }
    
    list(name=x[["name"]],
         union=m_union,
         intersection=m_intersect,
         vol_u=vol_u,
         vol_i=vol_i)    
}

get_mesh_ui <- function(x, boov=FALSE) {
    pairL <- get_mesh_pairs(x)
    Map(get_mesh_ui_pair, pairL, boov=boov)
}

get_mesh_metro_pair <- function(x, chop=TRUE, ...) {
    metro <- vcgMetro(x[["mesh_1"]][["mesh"]]$getMesh(),
                      x[["mesh_2"]][["mesh"]]$getMesh(),
                      ...)
    
    if(chop) {
        metro[["distances1"]]    <- NULL
        metro[["distances2"]]    <- NULL
        metro[["forward_hist"]]  <- NULL
        metro[["backward_hist"]] <- NULL
    }
    
    metro[["mesh_1"]] <- metro[["mesh1"]]
    metro[["mesh_2"]] <- metro[["mesh2"]]
    metro[["mesh1"]]  <- NULL
    metro[["mesh2"]]  <- NULL
    metro[["name"]]   <- x[["name"]]
    metro[["group"]]  <- x[["group"]]
    metro
}

get_mesh_metro <- function(x, chop=TRUE, ...) {
    pairL <- get_mesh_pairs(x)
    Map(get_mesh_metro_pair, pairL, chop=chop, ...)
}

## distance measures, union, intersection for each mesh pair
get_mesh_agree_pair <- function(x, metro, ui, boov=FALSE, do_ui=FALSE, chop=TRUE, ...) {
    ## distance-based measures
    if(missing(metro)) {
        metro <- get_mesh_metro_pair(x, chop=chop, ...)
    }
    
    DCOM        <- sqrt(sum((x[["mesh_2"]][["centroid"]] -
                             x[["mesh_1"]][["centroid"]])^2))
    HD_forward  <- metro[["ForwardSampling"]][["maxdist"]]
    HD_backward <- metro[["BackwardSampling"]][["maxdist"]]

    if(is.finite(HD_forward) && is.finite(HD_backward)) {
        HD_max <- max(c(HD_forward, HD_backward))
        HD_avg <- (HD_forward + HD_backward) / 2
    } else {
        HD_max <- NA_real_
        HD_avg <- NA_real_
    }
    
    ## average surface distance based on weighted average of sampled distances
    ## not on actual vertex distances as stored in distances1, distances2
    n1 <- metro[["ForwardSampling"]][["nsamples"]]
    n2 <- metro[["BackwardSampling"]][["nsamples"]]
    if((n1 > 0L) && (n2 > 0L)) {
        w1   <- n1 / (n1+n2)
        w2   <- n2 / (n1+n2)
        ASD  <-      w1* metro[["ForwardSampling"]][["meandist"]]   +
                     w2* metro[["BackwardSampling"]][["meandist"]]
        RMSD <- sqrt(w1*(metro[["ForwardSampling"]][["RMSdist"]]^2) +
                     w2*(metro[["BackwardSampling"]][["RMSdist"]]^2))
    } else {
        ASD  <- NA_real_
        RMSD <- NA_real_
    }
    
    ## volume-overlap-based measures
    ## check if union/intersection are supplied
    if(missing(ui) && do_ui) {
        ui <- get_mesh_ui_pair(x, boov=boov)
    }
    
    vol_1 <- x[["mesh_1"]][["volume"]]
    vol_2 <- x[["mesh_2"]][["volume"]]
    
    if(do_ui && !is.null(ui) && !is.null(ui[["union"]]) && !is.null(ui[["intersection"]])) {
        vol_u <- ui[["vol_u"]]
        vol_i <- ui[["vol_i"]]
        JSC   <-   vol_i / vol_u
        DSC   <- 2*vol_i / (vol_1 + vol_2)
        
        ## TODO TP, FP, TN, FN
    } else {
        vol_u <- NA_real_
        vol_i <- NA_real_
        JSC   <- NA_real_
        DSC   <- NA_real_
    }
    
    data.frame(mesh_1=x[["mesh_1"]][["name"]],
               mesh_2=x[["mesh_2"]][["name"]],
               group =x[["group"]],
               vol_1 =vol_1,
               vol_2 =vol_2,
               vol_u =vol_u,
               vol_i =vol_i,
               DCOM  =DCOM,
               HD_max=HD_max,
               HD_avg=HD_avg,
               ASD   =ASD,
               RMSD  =RMSD,
               JSC   =JSC,
               DSC   =DSC)
}

get_mesh_agree <- function(x, boov=FALSE, do_ui=FALSE, chop=TRUE, ...) {
    pairL  <- get_mesh_pairs(x)
    metroL <- Map(get_mesh_metro_pair, pairL, chop=chop, ...)
    uiL    <- if(do_ui) {
        Map(get_mesh_ui_pair, pairL)
    } else {
        list(NULL)
    }

    agreeL <- Map(get_mesh_agree_pair,
                  pairL,
                  metro=metroL,
                  ui   =uiL,
                  boov =boov,
                  do_ui=do_ui,
                  chop =chop)

    d <- do.call("rbind", agreeL)
    rownames(d) <- NULL
    d
}

get_mesh_agree_long <- function(x) {
    vars_varying <- c("DCOM",
                      "HD_max", "HD_avg", "ASD", "RMSD",
                      "vol_u", "vol_i",
                      "JSC", "DSC")
    
    vars_id <- names(x)[!(names(x) %in% vars_varying)]
    
    dL <- reshape(x,
                  direction="long",
                  idvar    =vars_id,
                  varying  =vars_varying,
                  v.names  ="observed",
                  timevar  ="metric")
    
    rownames(dL) <- NULL
    
    dL[["metric"]] <- factor(dL[["metric"]],
                             levels=seq_along(vars_varying),
                             labels=vars_varying)
    
    dL
}

get_mesh_agree_aggr <- function(x, na.rm=FALSE) {
    d_agreeL <- get_mesh_agree_long(x)
    d_agreeL[["observed_ln"]] <- log(d_agreeL[["observed"]])
    
    d_mean   <- aggregate(observed    ~ group + metric, FUN=mean,   data=d_agreeL, na.rm=na.rm)
    d_median <- aggregate(observed    ~ group + metric, FUN=median, data=d_agreeL, na.rm=na.rm)
    d_sd     <- aggregate(observed    ~ group + metric, FUN=sd,     data=d_agreeL, na.rm=na.rm)
    d_var    <- aggregate(observed    ~ group + metric, FUN=var,    data=d_agreeL, na.rm=na.rm)
    d_varlog <- aggregate(observed_ln ~ group + metric, FUN=var,    data=d_agreeL, na.rm=na.rm)
    
    d_aggr <- Reduce(function(x, y) { suppressWarnings(merge(x, y, by=c("group", "metric"))) },
                     list(d_mean, d_median, d_sd, d_var, d_varlog))
    
    names(d_aggr)       <- c("group", "metric", "Mean", "Median", "SD", "VAR", "VAR_log")
    d_aggr[["CV"]]      <- d_aggr[["SD"]] / d_aggr[["Mean"]]
    d_aggr[["CV_ln"]]   <- sqrt(exp(d_aggr[["VAR_log"]]) - 1)
    d_aggr[["VAR"]]     <- NULL
    d_aggr[["VAR_log"]] <- NULL
    
    d_aggr
}

get_mesh_agree_aggr_long <- function(x) {
    vars_varying <- c("Mean", "Median", "SD", "CV", "CV_ln")
    vars_id      <- names(x)[!(names(x) %in% vars_varying)]
    
    dL <- reshape(x,
                  direction="long",
                  idvar    =vars_id,
                  varying  =vars_varying,
                  v.names  ="observed",
                  timevar  ="statistic")
    
    rownames(dL) <- NULL
    
    dL[["statistic"]] <- factor(dL[["statistic"]],
                                levels=seq_along(vars_varying),
                                labels=vars_varying)

    dL
}

meshL_to_observerL <- function(x) {
    ll <- Map(function(i, name) {
        setNames(list(i), name)
    }, x, names(x))
    
    setNames(ll, sprintf("Observer_%.2d", seq_along(ll)))
}

mesh3dL_to_cgalMeshL <- function(x) {
    convert_mesh_one <- function(y) {
        if(inherits(y[["mesh"]], "mesh3d")) {
            y[["mesh"]] <- cgalMeshes::cgalMesh$new(y[["mesh"]])
        }
        
        y
    }
    
    convert_meshL <- function(z) {
        lapply(z, convert_mesh_one)
    }
    
    lapply(x, convert_meshL)
}
