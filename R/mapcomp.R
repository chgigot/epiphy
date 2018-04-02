#------------------------------------------------------------------------------#
#' Map Comparison procedure.
#'
#' \code{mapcomp} procedure. Spatial analyses of ecological count data: A density map comparaison approach
#'
#' the Hellinger distance.
#'
#' @param data A data frame with 3 columns.
#'
#' @references
#'
#' Lavigne C, Ricci B, Franck P, Senoussi R. 2010. Spatial analyses of
#' ecological count data: A density map comparison approach. Basic and Applied
#' Ecology. 11:734–742.
#'
#' @examples
#' my_res <- mapcomp(codling_moths, 4, 9, edge_correction = FALSE)
#' plot(my_res)
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
mapcomp.data.frame <- function(data, delta, h, nperm = 100,
                               edge_correction = TRUE, threads = 1) {

    # data structure:
    # - 1st and 2nd columns: x and y coordinates, respectively.
    # - 3rd column: observed disease intensity data.
    stopifnot(ncol(data) == 3)
    colnames(data) <- c("x", "y", "i")
    # ^ Needed for ggplot2 figures and to simplify the code below.

    # delta = mesh size of G (delta = delta_min here)
    # Define the mesh G
    grid_inter <- mesh_intersect(data, delta_min = delta)
    sites_homogene <- data
    sites_homogene[, "i"] <- 1
    #browser()
    phs   <- p_hscaled(as.matrix(grid_inter), as.matrix(data), h, edge_correction)
    qhs   <- p_hscaled(as.matrix(grid_inter), as.matrix(sites_homogene), h, edge_correction)
    test  <- delta / sqrt(2) * sqrt( sum( ( sqrt(phs) - sqrt(qhs) )^2 ) )
    res <- pbapply::pbsapply(seq_len(nperm), function(i) {
        new_data <- data
        new_data[, 3] <- sample(new_data[, 3])
        #sub_mapcomp(as.matrix(grid_inter), as.matrix(data),
        #            as.matrix(sites_homogene), delta, h, edge_correction)[["test"]]
        phs   <- p_hscaled(as.matrix(grid_inter), as.matrix(new_data), h, edge_correction)
        qhs   <- p_hscaled(as.matrix(grid_inter), as.matrix(sites_homogene), h, edge_correction)
        test  <- delta / sqrt(2) * sqrt( sum( ( sqrt(phs) - sqrt(qhs) )^2 ) )
        test
    }, cl = threads)
    coord <- data.frame(grid_inter, phs = phs)
    structure(list(data = data, coord = coord, test = test,
                   pval = (sum(res > test) + 1) / (nperm + 1)), # TODO: To check
              class = "mapcomp")
}

#------------------------------------------------------------------------------#
#' @rdname mapcomp
#' @export
#------------------------------------------------------------------------------#
mapcomp.matrix <- function(data, delta, h, nperm = 100,
                           edge_correction = TRUE, threads = 1) {
    mapcomp.data.frame(as.data.frame(data), delta, h, nperm, edge_correction,
                     threads)
}

#------------------------------------------------------------------------------#
#' @rdname mapcomp
#' @export
#------------------------------------------------------------------------------#
mapcomp.count <- function(data, delta, h, nperm = 100,
                          edge_correction = TRUE, threads = 1) {
    mapped_data <- map_data(data)
    mapped_data <- mapped_data[, c("x", "y", "i")] # no t
    mapcomp.data.frame(mapped_data, delta, h, nperm, edge_correction, threads)
}

#------------------------------------------------------------------------------#
#' @rdname mapcomp
#' @export
#------------------------------------------------------------------------------#
mapcomp.incidence <- function(data, delta, h, nperm = 100,
                              edge_correction = TRUE, threads = 1) {
    mapped_data <- map_data(data)
    mapped_data <- mapped_data[, c("x", "y", "i")] # no t, no n
    mapcomp.data.frame(mapped_data, delta, h, nperm, edge_correction, threads)
}


#==============================================================================#
# Print, summary and plot
#==============================================================================#

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
plot.mapcomp <- function(x, ...) {
    gg <- ggplot() +
        geom_raster(inherit.aes = FALSE, data = x$coord,
                    aes(x, y, fill = phs)) +
        geom_contour(inherit.aes = FALSE, data = x$coord,
                     aes(x, y, z = phs), size = 0.6, color = "black") +
        geom_point(inherit.aes = FALSE, data = x$data,
                   aes(x, y, size = i)) +
        scale_fill_gradient(low = "white", high = "red") +
        theme_bw()
    print(gg)
    invisible(NULL)
}


#==============================================================================#
# Utilities
#==============================================================================#

mesh_intersect <- function(sites, delta_min, delta_max = 2 * delta_min, ...,
                           threads = 1) {#parallel::detectCores()) {
    xrange <- range(sites[, "x"]) + c(-delta_min, delta_min)
    yrange <- range(sites[, "y"]) + c(-delta_min, delta_min)
    x <- seq(xrange[1], xrange[2], by = delta_min)
    y <- seq(yrange[1], yrange[2], by = delta_min)
    mesh_coords <- expand.grid(x = x, y = y, KEEP.OUT.ATTRS = FALSE)
    mesh_coords
    ###n <- nrow(sites)
    ###idx <- numeric(0)
    ###invisible(pbapply::pblapply(seq_len(nrow(mesh_coords)), function(i) {
    ###    out <- rowSums(abs(stack_rep(mesh_coords[i, ], n) - sites[, 1:2]) <
    ###                       stack_rep(c(delta_max, delta_max), n))
    ###    if (any(out == 2)) idx <<- c(idx, i) # Means if at least one site is in the neighboorhood
    ###}, cl = threads))
    ###mesh_coords[idx, ]
}


# Test, distance to homogeneity
# p: density of counts
# q: sampling density

#stack_rep <- function(x, n = 1) {
#    if (!is.vector(x)) x <- unlist(x) # C'est peutêtre une data.frame
#    ncol <- length(x)
#    matrix(rep(x, n), ncol = ncol, byrow = TRUE)
#}

