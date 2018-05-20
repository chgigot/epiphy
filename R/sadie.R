#------------------------------------------------------------------------------#
#' Spatial Analysis by Distance IndicEs (SADIE).
#'
#' \code{sadie} performs the SADIE procedure. It computes different indices and
#' probabilities based on the distance to regularity for the observed spatial
#' pattern and a specified number of random permutations of this pattern. Both
#' kind of clustering indices described by Perry et al. (1999) and Li et al.
#' (2012) can be computed.
#'
#' By convention in the SADIE procedure, clustering indices for a donor unit
#' (outflow) and a receiver unit (inflow) are positive and negative in sign,
#' respectively.
#'
#' @param data A data frame or a matrix with only three columns: the two first
#'     ones must be the x and y coordinates of the sampling units, and the last
#'     one, the corresponding disease intensity observations. It can also be a
#'     \code{\link{count}} or an \code{\link{incidence}} object.
#' @param index The index to be calculated: "Perry", "Li-Madden-Xu" or "all".
#'     By default, only Perry's index is computed for each sampling unit.
#' @param nperm Number of random permutations to assess probabilities.
#' @param seed Fixed seed to be used for randomizations (only useful for
#'     checking purposes). Not fixed by default (= NULL).
#' @param threads Number of threads to perform the computations.
#' @param method Method for the transportation algorithm.
#' @param verbose Explain what is being done (TRUE by default).
#' @param ... Additional arguments to be passed to other methods.
#'
#' @references
#'
#' Perry JN. 1995. Spatial analysis by distance indices. Journal of Animal
#' Ecology 64, 303–314. \href{http://dx.doi.org/10.2307/5892}{doi:10.2307/5892}
#'
#' Perry JN, Winder L, Holland JM, Alston RD. 1999. Red–blue plots for detecting
#' clusters in count data. Ecology Letters 2, 106–113.
#' \href{http://dx.doi.org/10.1046/j.1461-0248.1999.22057.x}{doi:10.1046/j.1461-0248.1999.22057.x}
#'
#' Li B, Madden LV, Xu X. 2012. Spatial analysis by distance indices: an
#' alternative local clustering index for studying spatial patterns. Methods in
#' Ecology and Evolution 3, 368–377.
#' \href{http://dx.doi.org/10.1111/j.2041-210X.2011.00165.x}{doi:10.1111/j.2041-210X.2011.00165.x}
#'
#' @examples
#' set.seed(123)
#' # Create an intensity object:
#' my_count <- count(aphids, mapping(x = xm, y = ym))
#' # Only compute Perry's indices:
#' my_res <- sadie(my_count)
#' my_res
#' summary(my_res)
#' plot(my_res)
#' plot(my_res, isoclines = TRUE)
#'
#' set.seed(123)
#' # Compute both Perry's and Li-Madden-Xu's indices (using multithreading):
#' my_res <- sadie(my_count, index = "all", threads = 2, nperm = 20)
#' my_res
#' summary(my_res)
#' plot(my_res) # Identical to: plot(my_res, index = "Perry")
#' plot(my_res, index = "Li-Madden-Xu")
#'
#' set.seed(123)
#' # Using usual data frames instead of intensity objects:
#' my_df <- aphids[, c("xm", "ym", "i")]
#' sadie(my_df)
#'
#' @name sadie
#' @export
#------------------------------------------------------------------------------#
sadie <- function(data, ...) UseMethod("sadie")

#------------------------------------------------------------------------------#
#' @rdname sadie
#' @method sadie data.frame
#' @export
#------------------------------------------------------------------------------#
sadie.data.frame <- function(data, index = c("Perry", "Li-Madden-Xu", "all"),
                             nperm = 100, seed = NULL, threads = 1, ...,
                             method = "shortsimplex", verbose = TRUE) {

    dots <- list(...)
    if (is.null(call <- dots[["call"]])) {
        call <- match.call()
    }

    index <- match.arg(index)
    # data structure:
    # - 1st and 2nd columns: x and y coordinates, respectively.
    # - 3rd column: observed disease intensity data.
    stopifnot(ncol(data) == 3)
    colnames(data) <- c("x", "y", "i")
    # ^ Needed for ggplot2 figures and to simplify the code below.
    cost <- as.matrix(dist(data[c("x", "y")]))

    # More interesting
    if (!is.null(seed)) set.seed(seed)

    N <- nrow(data)
    data[["i"]] <- as.numeric(data[["i"]]) # Just in case it's an integer, to have the same nature for data and flat_data

    # Create a homogeneous (flatten) version of data
    flat_data <- rep(mean(data[["i"]]), N)

    # Use optimal transportation algorithm
    opt_transport <- wrap_transport(data[["i"]], flat_data, cost, method)
    opt_transport <- as.matrix(opt_transport, N)

    # Computation the costs
    # Due to the convention only positive values (donor units) are added up
    # to compute Da. Da = total observed distance (or total cost).
    cost_flows <- costTotCPP(opt_transport, cost)
    Da <- sum(cost_flows[cost_flows >= 0]) # Donors

    # Deal with indices
    if (index == "all") index <- c("Perry", "Li-Madden-Xu")
    info_P   <- list(clust = NA_real_, Ea = NA_real_)
    info_LMX <- list(clust = NA_real_, Dis_all = matrix(NA_real_))
    idx_P    <- NA_real_
    idx_LMX  <- NA_real_
    Ea       <- NA_real_ # Ea = randomized distances
    Ia       <- NA_real_ # Da = observed,
    Pa       <- NA_real_
    new_prob <- NA_real_
    Dis_all  <- NA_real_
    if (any(index == "Perry")) {
        info_P <- clust_P(opt_transport, cost, N, nperm,
                          start = data[["i"]], end = flat_data,
                          method, cost_flows, threads, verbose)
        idx_P <- info_P[["clust"]]
        Ea    <- info_P[["Ea"]]
        Ia    <- Da / mean(Ea)
        Pa    <- sum(Ea > Da) / length(Ea)
    } else {
        warning("Ia and Pa cannot be computed if Perry's indices are not.")
    }
    if (any(index == "Li-Madden-Xu")) {
        info_LMX <- clust_LMX(opt_transport, cost, N, nperm,
                              start = data[["i"]], end = flat_data,
                              method, cost_flows, threads, verbose)
        idx_LMX  <- info_LMX[["clust"]]
        new_prob <- info_LMX[["prob"]]
        Dis_all  <- info_LMX[["Dis_all"]]
    }

    info_clust <- data.frame(data,
                             cost_flows = cost_flows,
                             idx_P      = idx_P,
                             idx_LMX    = idx_LMX,
                             prob       = new_prob)

    # Summary indices
    summary_idx <- data.frame(overall = c(mean(abs(idx_P)),
                                          mean(abs(idx_LMX))),
                              inflow  = c(mean(idx_P[idx_P < 0]),
                                          mean(idx_LMX[idx_LMX < 0])),
                              outflow = c(mean(idx_P[idx_P > 0]),
                                          mean(idx_LMX[idx_LMX > 0])))
    rownames(summary_idx) <- c("Perry's index", "Li-Madden-Xu's index")

    # Returns:
    res <- list(info_clust = info_clust,
                Da = Da, # Da = total observed distance (or total cost).
                Ia = Ia,
                Pa = Pa,
                Ea = Ea,
                summary_idx = summary_idx,
                nperm = nperm,
                seed  = seed)
    attr(res, "class") <- "sadie"
    attr(res, "call")  <- call
    res
}

#------------------------------------------------------------------------------#
#' @rdname sadie
#' @export
#------------------------------------------------------------------------------#
sadie.matrix <- function(data, index = c("Perry", "Li-Madden-Xu", "all"),
                         nperm = 100, seed = NULL, threads = 1, ...,
                         method = "shortsimplex", verbose = TRUE) {
    sadie.data.frame(as.data.frame(data), index, nperm, seed, threads, ...,
                     method = method, verbose = verbose, call = match.call())
}

#------------------------------------------------------------------------------#
#' @rdname sadie
#' @export
#------------------------------------------------------------------------------#
sadie.count <- function(data, index = c("Perry", "Li-Madden-Xu", "all"),
                        nperm = 100, seed = NULL, threads = 1, ...,
                        method = "shortsimplex", verbose = TRUE) {
    mapped_data <- map_data(data)
    mapped_data <- mapped_data[, c("x", "y", "i")] # no t
    sadie.data.frame(mapped_data, index, nperm, seed, threads, ...,
                     method = method, verbose = verbose, call = match.call())
}

#------------------------------------------------------------------------------#
#' @rdname sadie
#' @export
#------------------------------------------------------------------------------#
sadie.incidence <- function(data, index = c("Perry", "Li-Madden-Xu", "all"),
                        nperm = 100, seed = NULL, threads = 1, ...,
                        method = "shortsimplex", verbose = TRUE) {
    mapped_data <- map_data(data)
    mapped_data <- mapped_data[, c("x", "y", "i")] # no t, no n
    sadie.data.frame(mapped_data, index, nperm, seed, threads, ...,
                     method = method, verbose = verbose, call = match.call())
}


#==============================================================================#
# Print, summary and plot
#==============================================================================#

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.sadie <- function(x, ...) {
    cat("Spatial Analysis by Distance IndicEs (sadie)\n")
    cat("\nCall:\n")
    print(attr(x, "call"))
    cat("\nIa: ", format(x$Ia, digits = 1, nsmall = 4),
        " (Pa = ", format.pval(x$Pa), ")\n\n", sep = "")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
summary.sadie <- function(object, ...) {
    res <- object
    class(res) <- "summary.sadie"
    res
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.summary.sadie <- function(x, ...) {
    cat("\nCall:\n")
    print(attr(x, "call"))
    n <- 6L
    cat("\nFirst ", n, " rows of clustering indices:\n", sep = "") # Cf Image pkg for a nice display
    print(head(x$info_clust, n = n))
    #printCoefmat(x$info_clust)
    #cat("number of permutations: ", object$nperm, "\n",
    #    "random seed: ", object$rseed, "\n", sep = "")
    cat("\nSummary indices:\n")
    print(x$summary_idx)
    cat("\nMain outputs:")
    cat("\nIa: ", format(x$Ia, digits = 1, nsmall = 4),
        " (Pa = ", format.pval(x$Pa), ")\n", sep = "")
    cat("\n'Total cost': ", x$Da, "\n",
        "Number of permutations: ", x$nperm, "\n", sep = "")
    if (!is.null(x$seed)) { # If fixed seed
        cat("Fixed seed: ", x$seed, "\n", sep = "")
    }
    cat("\n")
}

#------------------------------------------------------------------------------#
# TODO: Doc about thresholds, etc. isoclines good when well sampled, the grid!
#' @export
#------------------------------------------------------------------------------#
plot.sadie <- function(x, ..., index = c("Perry", "Li-Madden-Xu"),
                       isoclines = FALSE, resolution = 100, bins = 5,
                       thresholds = c(-1.5, 1.5), conf.level = 0.95,
                       point_size = c("radius", "area")) {

    index      <- match.arg(index)
    point_size <- match.arg(point_size)
    data_clust <- x$info_clust # To make it simplier.

    idx <- switch(index, "Perry" = "idx_P", "Li-Madden-Xu" = "idx_LMX")
    data_clust$col <- switch(index,
        "Perry" = {
            vapply(data_clust[[idx]], function(x) {
                if      (x <  thresholds[1L]) "blue"  # Low values
                else if (x <= thresholds[2L]) "white" # Intermediate values
                else                          "red"   # High values
            }, character(1L))
        },
        "Li-Madden-Xu" = {
            vapply(seq_len(nrow(data_clust)), function(i1) {
                if (data_clust[["prob"]][i1] <= (1 - conf.level)) {
                    if      (data_clust[[idx]][i1] <  0) "blue" # Low values
                    else if (data_clust[[idx]][i1] >= 0) "red"  # High values
                } else {
                    "white" # Intermediate values
                }
            }, character(1L))
        })

    gg <- ggplot()
    if (isoclines) {
        # Prepare interpolated data:
        data_loess <- stats::loess(as.formula(paste0(idx, " ~ x * y")),
                                   data = data_clust, degree = 2, span = 0.2)
        input_val <- expand.grid(
            x = seq(min(data_clust$x), max(data_clust$x), length = resolution),
            y = seq(min(data_clust$y), max(data_clust$y), length = resolution))
        interpolated   <- predict(data_loess, input_val)
        data_landscape <- data.frame(input_val, z = as.vector(interpolated))

        # Build ggplot figure:
        gg <- gg + geom_raster(data = data_landscape, aes(x, y, fill = z))
        gg <- gg + geom_contour(data = data_landscape, aes(x, y, z = z),
                                bins = bins, size = 0.6, color = "black")
        # Note that it is not possible in current ggplot2 implementation to use
        # two different colour scales in the same figure.
        gg <- gg + geom_point(data = data_clust,
                              aes_(quote(x), quote(y),
                                   size = call("abs", as.name(idx))),
                              colour = "black", fill = "black", pch = 21)
        switch(point_size,
            "area" = {
                gg <- gg + scale_size("Absolute\nindex", range = c(0, 10))
            },
            "radius" = {
                gg <- gg + scale_radius("Absolute\nindex", range = c(0, 10))
            }
        )
        gg <- gg + scale_fill_gradientn("Interpolated\nindex",
                                        colours = terrain.colors(10))
    } else {
        # Build ggplot figure:
        gg <- gg + geom_point(data = data_clust,
                              aes_(quote(x), quote(y), fill = quote(col),
                                   size = call("abs", as.name(idx))),
                              colour = "black",  pch = 21)
        switch(point_size,
            "area" = {
                gg <- gg + scale_size("Absolute\nindex", range = c(0, 10))
             },
             "radius" = {
                gg <- gg + scale_radius("Absolute\nindex", range = c(0, 10))
            }
        )
        switch(index,
            "Perry" = {
                gg <- gg + scale_fill_manual("Index\nthreshold",
                    breaks = c("red", "white", "blue"),
                    labels = c("red"   = paste0("> ", thresholds[2L]),
                               "white" = paste0("Between ", thresholds[1L],
                                                "\nand ", thresholds[2L]),
                               "blue" = paste0("< ", thresholds[1L])),
                    values = c("red" = "red",
                               "white" = "white",
                               "blue" = "blue"))
            },
            "Li-Madden-Xu" = {
                gg <- gg + scale_fill_manual("Index\nthreshold",
                    breaks = c("red", "white", "blue"),
                    labels = c("red"   = "Sign. > 0",
                               "white" = "No sign.\ndifferent\nfrom 0",
                               "blue" = "Sign. < 0"),
                    values = c("red" = "red",
                               "white" = "white",
                               "blue" = "blue"))
            })
    }
    gg <- gg + theme_bw()
    print(gg)
    invisible(NULL)
}


#==============================================================================#
# Utilities
#==============================================================================#

#------------------------------------------------------------------------------#
# Wrapper for the transport::transport function.
#------------------------------------------------------------------------------#
wrap_transport <- function(start, end, cost, method = "shortsimplex",
                           control = list(), ...) {
    res <- transport::transport(start, end, costm = cost, method = method,
                                control = control, ...)
    class(res) <- c("transport", class(res))
    res
}

#------------------------------------------------------------------------------#
#' @method as.matrix transport
#------------------------------------------------------------------------------#
as.matrix.transport <- function(x, dim_mat) {
    as_matrix_transport(x, dim_mat)
}

#------------------------------------------------------------------------------#
# Computation of clustering indices from Perry et al. (1999).
#------------------------------------------------------------------------------#
clust_P <- function(flow, cost, dim_mat, nperm, start, end,
                    method = "shortsimplex", cost_flows, threads,
                    verbose = TRUE) {
    if (!verbose) {
        op <- pbapply::pboptions(type = "none")
    } else {
        cat("Computation of Perry's indices:\n")
    }

    # Randomization procedure:
    idx <- seq_len(nrow(flow)) # or: seq_len(length(start))
    randomizations <- pbapply::pblapply(seq_len(nperm), function(i1) {
        rand_idx       <- sample(idx)
        new_start      <- start[rand_idx]
        opt_transport  <- wrap_transport(new_start, end, cost, method)
        opt_transport  <- as.matrix(opt_transport, dim_mat)
        idx_same_count <- lapply(idx, function(i2) which(new_start == start[i2]))
        res_Yci <- vapply(seq_len(length(idx_same_count)), function(i2) {
            mean(vapply(seq_len(length(idx_same_count[[i2]])),
                function(x) {
                    costTotiCPP(idx_same_count[[i2]][x], opt_transport, cost)
                }, numeric(1L)))
        }, numeric(1L))
        res_Yii <- vapply(idx,
            function(x) {
                costTotiCPP(x, opt_transport, cost)
            }, numeric(1L))
        # Due to the convention, only positive values (donor units) are added up
        # to compute Ea.
        cost_flows_perm <- costTotCPP(opt_transport, cost)
        Ea <- sum(cost_flows_perm[cost_flows_perm >= 0]) # Donors
        list(flow = opt_transport,
             idx  = rand_idx,
             Yci  = res_Yci,
             Yii  = res_Yii,
             cost_flows_perm = Ea)
    }, cl = threads)

    # Compute indices:
    Ycs <- simplify2array(lapply(randomizations, function(x) x[["Yci"]]))
    Ycs <- rowMeans(abs(Ycs))
    Yis <- simplify2array(lapply(randomizations, function(x) x[["Yii"]]))
    Yis <- rowMeans(abs(Yis))
    Y0  <- mean(Ycs) # Must be equal to: mean(Yis)
    Yi  <- vapply(idx,
        function(x) {
            costTotiCPP(x, flow, cost)
        }, numeric(1L))
    clust <- ((Yi * Y0) / (Yis * Ycs))
    Ea <- unlist(lapply(randomizations, function(x) x[["cost_flows_perm"]]))

    if (!verbose) {
        pbapply::pboptions(op)
    }

    # Returns:
    list(clust = clust,
         Ea    = Ea)
}

#------------------------------------------------------------------------------#
# Computation of clustering indices from Li, Madden and Xu (2012).
#------------------------------------------------------------------------------#
clust_LMX <- function(flow, cost, dim_mat, nperm, start, end,
                      method = "shortsimplex", cost_flows, threads,
                      verbose = TRUE) {

    if (!verbose) {
        op <- pbapply::pboptions(type = "none")
    } else {
        cat("Computation of Li-Madden-Xu's indices:\n")
    }

    # Randomization procedure:
    idx <- seq_len(nrow(flow)) # or: seq_len(length(start))
    randomizations <- pbapply::pblapply(seq_len(nperm), function(i1) {
        vapply(idx, function(i2) {
            sub_idx       <- idx[idx != i2]
            rand_idx      <- sample(sub_idx)
            rand_idx      <- append(rand_idx, i2, after = i2 - 1)
            new_start     <- start[rand_idx]
            opt_transport <- wrap_transport(new_start, end, cost, method)
            opt_transport <- as.matrix(opt_transport, dim_mat)
            costTotiCPP(i2, opt_transport, cost)
        }, numeric(1L))
    }, cl = threads)

    # Compute indices:
    randomizations <- simplify2array(randomizations)
    Dis_bar        <- rowMeans(cbind(randomizations, cost_flows))
    clust   <- cost_flows / abs(Dis_bar)
    clust_i <- randomizations / abs(Dis_bar)
    prob    <- vapply(seq_len(length(clust)), function(i1) {
        if (clust[i1] >= 0) rpos <- (clust_i[i1, ] > clust[i1]) # Donor
        else                rpos <- (clust_i[i1, ] < clust[i1]) # Receiver
        rpos <- sum(rpos) + 1
        # Trick above:
        # (1) sum(T, T, F) gives 2.
        # (2) +1 b/c clust is +1 away from the nearest lower or upper clust_i.
        rpos / (nperm + 1)
    }, numeric(1L))

    if (!verbose) {
        pbapply::pboptions(op)
    }

    # Returns:
    list(clust   = clust,
         prob    = prob)
}







