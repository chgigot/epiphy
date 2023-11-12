#------------------------------------------------------------------------------#
#' Map Comparison procedure.
#'
#' \code{mapcomp} performs a spatial pattern analysis based on the calculation
#' of a formal distance (the Hellinger distance) between the density map of
#' count or incidence data, and the density map of sampling effort. Statistical
#' tests of spatial homogeneity are based on permutations across sampling sites
#' and on valuable properties of the Hellinger distance.
#'
#' @param data A data frame or a matrix with only three columns: the two first
#'     ones must be the x and y coordinates of the sampling units, and the last
#'     one, the corresponding disease intensity observations. It can also be a
#'     \code{\link{count}} or an \code{\link{incidence}} object.
#'
#' @param bandwidth Bandwidth parameter for smoothing. It allows to test the
#'     spatial extent of heterogeneity if any.
#' @param delta Mesh size of the grid over the geographical domain of the
#'     sampling units used to compute the integral Hellinger distance between
#'     the probability density function of observations and the probability
#'     density function of sampling effort.
#' @param edge_correction Apply edge correction to account for the fact that
#'     bordering points intrinsically suffer from a lack of neighboring
#'     observation sites. FALSE by default.
#' @param nperm Number of random permutations to assess probabilities.
#' @param threads Number of threads to perform the computations.
#' @param verbose Explain what is being done (TRUE by default).
#' @param ... Additional arguments to be passed to other methods.
#'
#' @returns
#' An object of class \code{mapcomp}, which is a list containing the following
#' components:
#' \tabular{ll}{
#'     \code{data}       \tab The input data. \cr
#'     \code{coord}      \tab The coordinates and normalized intensity for each point of the full grid. \cr
#'     \code{object}     \tab The class of \code{data}. \cr
#'     \code{bandwidth}  \tab The \code{bandwidth} parameter. \cr
#'     \code{stat, pval} \tab The statistic and corresponding p-value (see references for more details). \cr
#' }
#'
#' @references
#'
#' Lavigne C, Ricci B, Franck P, Senoussi R. 2010. Spatial analyses of
#' ecological count data: A density map comparison approach. Basic and Applied
#' Ecology. 11:734–742.
#'
#' @examples
#' set.seed(123)
#' my_res <- mapcomp(codling_moths, delta = 1, bandwidth = 11,
#'                   edge_correction = FALSE, nperm = 20)
#' my_res
#' plot(my_res)
#'
#' set.seed(123)
#' my_count <- count(codling_moths, mapping(x = xm, y = ym))
#' my_res <- mapcomp(my_count, delta = 1, bandwidth = 11,
#'                   edge_correction = FALSE, nperm = 20)
#' my_res
#' plot(my_res, bins = 10)
#'
#' @name mapcomp
#' @export
#------------------------------------------------------------------------------#
mapcomp <- function(data, ...) UseMethod("mapcomp")

#------------------------------------------------------------------------------#
#' @rdname mapcomp
#' @method mapcomp data.frame
#' @export
#------------------------------------------------------------------------------#
mapcomp.data.frame <- function(data, delta, bandwidth, nperm = 100,
                               edge_correction = FALSE, threads = 1,
                               verbose = TRUE, ...) {

    if (!verbose) {
        op <- pbapply::pboptions(type = "none")
    }

    dots <- list(...)
    if (is.null(call <- dots[["call"]])) {
        call <- match.call()
    }
    # data structure:
    # - 1st and 2nd columns: x and y coordinates, respectively.
    # - 3rd column: observed disease intensity data.
    stopifnot(ncol(data) == 3)
    colnames(data) <- c("x", "y", "i")
    # ^ Needed for ggplot2 figures and to simplify the code below.

    # delta = mesh size of G (delta = delta_min here)
    # Define the mesh G
    grid_inter <- mesh_intersect(data, delta_min = delta)
    flat_data  <- data
    flat_data[, "i"] <- 1

    sub_mapcomp <- function(data, flat_data, grid_inter, delta, bandwidth,
                            edge_correction) {
        phs  <- p_hscaled(grid_inter, data, bandwidth, edge_correction)
        qhs  <- p_hscaled(grid_inter, flat_data, bandwidth, edge_correction)
        stat <- delta / sqrt(2) * sqrt(sum((sqrt(phs) - sqrt(qhs))^2))
        list(phs = phs, qhs = qhs, stat = stat)
    }

    res <- sub_mapcomp(data, flat_data, grid_inter, delta, bandwidth,
                       edge_correction)

    randomizations <- pbapply::pbsapply(seq_len(nperm), function(i) {
        new_data <- data
        new_data[, 3] <- sample(new_data[, 3])
        res <- sub_mapcomp(new_data, flat_data, grid_inter,
                    delta, bandwidth, edge_correction)
        res[["stat"]]
    }, cl = threads)

    coord <- data.frame(grid_inter, phs = res[["phs"]])
    res <- list(object = class(data),
                bandwidth = bandwidth,
                data  = data,
                coord = coord,
                stat  = res[["stat"]],
                pval  = (sum(randomizations > res[["stat"]]) + 1) / (nperm + 1)) # To double check
    attr(res, "class") <- "mapcomp"
    attr(res, "call")  <- call

    if (!verbose) {
        pbapply::pboptions(op)
    }

    res
}

#------------------------------------------------------------------------------#
#' @rdname mapcomp
#' @export
#------------------------------------------------------------------------------#
mapcomp.matrix <- function(data, delta, bandwidth, nperm = 100,
                           edge_correction = FALSE, threads = 1,
                           verbose = TRUE, ...) {
    mapcomp.data.frame(as.data.frame(data), delta, bandwidth, nperm,
                       edge_correction, threads, verbose, ...,
                       call = match.call())
}

#------------------------------------------------------------------------------#
#' @rdname mapcomp
#' @export
#------------------------------------------------------------------------------#
mapcomp.count <- function(data, delta, bandwidth, nperm = 100,
                          edge_correction = FALSE, threads = 1, verbose = TRUE,
                          ...) {
    mapped_data <- map_data(data)
    mapped_data <- mapped_data[, c("x", "y", "i")] # no t
    mapcomp.data.frame(mapped_data, delta, bandwidth, nperm, edge_correction,
                       threads, verbose, ..., call = match.call())
}

#------------------------------------------------------------------------------#
#' @rdname mapcomp
#' @export
#------------------------------------------------------------------------------#
mapcomp.incidence <- function(data, delta, bandwidth, nperm = 100,
                              edge_correction = FALSE, threads = 1,
                              verbose = TRUE, ...) {
    mapped_data <- map_data(data)
    mapped_data <- mapped_data[, c("x", "y", "i")] # no t, no n
    mapcomp.data.frame(mapped_data, delta, bandwidth, nperm, edge_correction,
                       threads, verbose, ..., call = match.call())
}


#==============================================================================#
# Print, summary and plot
#==============================================================================#

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.mapcomp <- function(x, ...) {
    cat("Map Comparison analysis (mapcomp)\n")
    cat("\nCall:\n")
    print(attr(x, "call"))
    cat("\nStat: ", format(x[["stat"]], digits = 1, nsmall = 4),
        " (P = ", format.pval(x[["pval"]]), ")\n\n", sep = "")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
plot.mapcomp <- function(x, bins = 5,...) {
    gg <- ggplot()
    gg <- gg + geom_raster(inherit.aes = FALSE, data = x$coord,
                           aes(x, y, fill = phs))
    gg <- gg + geom_contour(inherit.aes = FALSE, data = x$coord,
                            aes(x, y, z = phs),
                            bins = bins, size = 0.6, color = "black")
    gg <- gg + geom_point(inherit.aes = FALSE, data = x$data,
                          aes(x, y, size = i))
    gg <- gg + scale_size_continuous("Observed\nintensity")
    gg <- gg + scale_fill_gradient(paste0("Theoretical\nnormalized\n",
                                          "intensity for a\n",
                                          "bandwidth\nh = ", x[["bandwidth"]]),
                                   low = "white", high = "red")
    gg <- gg + theme_bw()
    print(gg)
    invisible(NULL)
}


#==============================================================================#
# Utilities
#==============================================================================#

mesh_intersect <- function(sites, delta_min, delta_max = 2 * delta_min, ...,
                           threads = 1) {
    xrange <- range(sites[, "x"]) + c(-delta_min, delta_min)
    yrange <- range(sites[, "y"]) + c(-delta_min, delta_min)
    expand.grid(x = seq(xrange[1], xrange[2], by = delta_min),
                y = seq(yrange[1], yrange[2], by = delta_min),
                KEEP.OUT.ATTRS = FALSE)
}



