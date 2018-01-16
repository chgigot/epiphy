#------------------------------------------------------------------------------#
#' @include utils.R
#' @include intensity-classes.R
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' Spatial hierarchy analysis
#'
#' TODO
#'
#' The manner in which the data are collected provides information about
#' aggregation of disease at different levels in a spatial hierarchy (Hughes et
#' al. 1997). For example, a sampling unit (upper level) can be reported as
#' "healthy", if no diseased leaves (lower level) were found within the sampling
#' unit.
#'
#' In a pairwise comparison between levels, the probability that an individual
#' at the lower hierarchical level is diseased is denoted plow, and phigh refers
#' to the probability of disease at the higher level. The relationship between
#' these two probabilities can be written as
#'
#' phigh = 1 - (1 - plow)^nu
#'
#' where n is a parameter ranging from 0 to the corresponding number of
#' individuals at the hierarchical level referenced by plow. If the value of n
#' is equal to the number of individuals at the lower hierarchical level
#' contained in a unit of the higher level (n low ), this suggests that there is
#' no aggregation of disease incidence at the lower level. Conversely, a value
#' of n less than n low is indicative of aggregation at that level. The value of
#' n can be interpreted as an effective sample size (Hughes and Gottwald 1999;
#' Madden and Hughes 1999) in the statistical sense that its value indicates the
#' number of independent pieces of information at the lower level. Here, the
#' effective sample size concerns the equating of the zero-term of the binomial
#' distribution with the zero-term of an overdispersed distribution, as
#' described in Madden and Hughes (1999). Using the complementary log-log
#' transformation, CLL(x) = ln(-ln(1-x)), one can rewrite the Equation 5 as
#' follows (Madden et al. 2007):
#'
#' CLL(phigh) = ln(nu) + CLL(plow)
#'
#' from which the value of ln(n) can be obtained as the intercept of a linear
#' regression when the slope is constrained to 1.
#'
#' @param low An list of \code{intensity} objects.
#' @param high An list of \code{intensity} objects.
#'
#' @examples
#' my_data_low  <- incidence(tomato_tswv$field_1928)
#' # TODO: 2 bugs to correct here (before)
#' # my_data_low  <- split(my_data_low, by = "t")[[1]]
#' # my_data_high <- clump(my_data_low, unit_size = c(x = 3, y = 3))
#'
#' my_data_low <- incidence(tomato_tswv$field_1929)
#' my_data_low <- clump(my_data_low, c(x = 3, y = 3))
#' my_data_high <- my_data_low
#' my_data_high$data$n <- 1
#' my_data_high$data$i <- ifelse(my_data_high$data$i > 0, 1, 0)
#' my_data_low  <- split(my_data_low, by = "t")
#' my_data_high <- split(my_data_high, by = "t")
#' res <- spatial_hier(my_data_low, my_data_high)
#'
#' res
#' summary(res)
#' plot(res)
#'
#' @export
#------------------------------------------------------------------------------#
spatial_hier <- function(low, high) {

    call   <- match.call()
    # Checks and variable allocation:
    if (missing(low) || missing(high)) {
        stop("Both 'low' and 'high' must be provided.")
    }
    object_class_low  <- unique(vapply(low, function(x) class(x)[[1L]],
                                       character(1L)))
    object_class_high <- unique(vapply(high, function(x) class(x)[[1L]],
                                       character(1L)))
    stopifnot(length(object_class_low) == 1 && length(object_class_high) == 1)
    stopifnot(object_class_low == object_class_high)
    if (length(low) != length(high)) {
        stop("'low' and 'high' lengths differ.")
    }

    # Get formatted data for low and high:
    data_low  <- get_fmt_obs(low, type = "incidence")
    data_high <- get_fmt_obs(high, type = "incidence")

    # Compute mean n:
    n_low  <- mean(vapply(data_low,  function(obj) mean(obj$n), numeric(1L)))
    n_high <- mean(vapply(data_high, function(obj) mean(obj$n), numeric(1L)))
    if (n_high != 1) {
        stop(paste0("Number of individuals per sampling units for 'high' must ",
                    "be equal to 1."))
    }

    # Compute mean p for each data set:
    p_low  <- vapply(data_low,  function(obj) mean(obj$p), numeric(1L))
    p_high <- vapply(data_high, function(obj) mean(obj$p), numeric(1L))

    # Perform the analysis
    coord_obs <- data.frame(x = p_low, y = p_high)
    data      <- cloglog(coord_obs) # TODO: Offer the possibility to use other log base
    data      <- data[with(data, is.finite(x) & is.finite(y)), ]
    if ((nr <- nrow(data)) < 1) {
        stop("Too few finite cloglog data to perform a regression.")
    } else if (nr < 5) {
        warning(paste0("Only ", nr, " cloglog data points were used to ",
                       "perform the regression."))
    }
    model     <- lm(y ~ offset(x), data = data)
    # ^ is the same as lm((y - x) ~ 1, ...), i.e. we just look for an intercept.
    nu        <- unname(exp(coef(model)))
    coord_the <- data.frame(x = p_low, y = 1 - (1 - p_low)^nu)

    # Return a spatial_hier object
    structure(list(call      = call, # TODO: Add more information about the transformation?
                   model     = model,
                   nu        = nu,    # TODO: create a param element like in power_law?
                   n         = n_low, # TODO: Where to put n???
                   coord_obs = coord_obs,
                   coord_the = coord_the),
              class = "spatial_hier")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.spatial_hier <- function(x, ...) { # TODO
    cat("# Spatial hierarchy analysis:\n")
    cat("nu = ", x$nu, "\n", sep = "")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
summary.spatial_hier <- function(object, ...) {
    # TODO: Bien regarder la structure d'un summary.lm pour bien comprendre
    # tous ses éléments.
    summary_model <- summary(object$model)
    summary_model$call <- object$call

    # TODO: Below, move to main function (like in power_law)??
    param   <- coef(summary(object$model))
    baseLog <- exp(1)
    nu      <- estimate_param(object$model, bquote(.(baseLog)^x1))
    param   <- rbind(param, unlist(nu))
    rownames(param) <- c("log_base(nu)", "nu")
    summary_model$coefficients <- param
    # ---

    structure(summary_model, class = "summary.spatial_hier")
}

#------------------------------------------------------------------------------#
#' @method print summary.spatial_hier
#' @export
#------------------------------------------------------------------------------#
print.summary.spatial_hier <- function(x, ...) stats:::print.summary.lm(x, ...)

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
plot.spatial_hier <- function(x, ..., scale = c("linear", "cloglog"),
                              observed = TRUE, model = TRUE, random = TRUE) {

    scale <- match.arg(scale)

    switch (scale,
        "linear" = {
            gg <- list(
                labs(x = expression(p[low]), y = expression(p[high])),
                scale_x_continuous(limits = c(0, 1)),
                scale_y_continuous(limits = c(0, 1)),
                switch(observed, geom_point(data = x$coord_obs, aes(x, y), ...)),
                switch(model, {
                    p_low  <- seq(0, 1, by = 0.01)
                    p_high <- 1 - (1 - p_low)^(x$nu)
                    geom_line(data = data.frame(x = p_low, y = p_high),
                              aes(x, y), ...)
                }),
                switch(random, {
                    p_low  <- seq(0, 1, by = 0.01)
                    p_high <- 1 - (1 - p_low)^(x$n)
                    geom_line(data = data.frame(x = p_low, y = p_high),
                              aes(x, y), linetype = "dashed", ...)
                }),
                theme_bw()
            )
        },
        "cloglog" = {
            gg <- list(
                labs(x = expression("cloglog(" * p[low] * ")"),
                     y = expression("cloglog(" * p[high] * ")")),
                switch(observed, geom_point(data = cloglog(x$coord_obs),
                                            aes(x, y), ...)),
                switch(model, {
                    p_low  <- seq(0, 1, by = 0.01)
                    p_high <- 1 - (1 - p_low)^(x$nu)
                    geom_line(data = cloglog(data.frame(x = p_low, y = p_high)),
                              aes(x, y), ...)
                }),
                switch(random, {
                    p_low  <- seq(0, 1, by = 0.001)
                    p_high <- 1 - (1 - p_low)^(x$n)
                    geom_line(data = cloglog(data.frame(x = p_low, y = p_high)),
                              aes(x, y), linetype = "dashed", ...)
                }),
                theme_bw()
            )
        }
    )
    ggplot() + gg
}

