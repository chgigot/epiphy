#------------------------------------------------------------------------------#
#' @include utils.R
#' @include intensity-classes.R
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' Taylor's and binary power laws.
#'
#' Assesses the overall degree of heterogeneity in a collection of data sets at
#' the sampling-unit scale.
#'
#' The power law describes the relationship between the observed variance of
#' individuals within a data set (\code{s^2}) and the corresponding variance
#' under the assumption of no aggregation (\code{s\'^2}). It can be expressed
#' under its logarithmic form as: \code{log(s^2) = log(a) + b log(Y)}, with:
#' \itemize{
#'     \item \code{Y = p} in the case of count data (Taylor's power law).
#'     \item \code{Y = p(1 - p)} in the case of incidence data (binary power law).
#' }
#' \code{p} corresponds to the mean proportion of recorded individuals in case
#' of incidence data, and the absolute value in case of count data.
#'
#' @param data A list of \code{intensity} objects (\code{count} or
#'     \code{incidence} objects).
#' @param log_base Logarithm base to be used.
#' @param ... Additional arguments to be passed to other methods.
#'
#' @return A \code{power_law} object.
#'
#' @examples
#' require(magrittr)
#' my_data <- do.call(c, lapply(citrus_ctv, function(citrus_field) {
#'    incidence(citrus_field) %>%
#'        clump(unit_size = c(x = 3, y = 3)) %>%
#'        split(by = "t")
#' }))
#' # my_data is a list of incidence object, each one corresponding to a given
#' # time at a given location.
#' my_power_law <- power_law(my_data)
#' my_power_law
#' summary(my_power_law)
#' plot(my_power_law) # Same as: plot(my_power_law, scale = "log")
#' plot(my_power_law, scale = "lin")
#'
#' @references
#'
#' Taylor LR. 1961. Aggregation, variance and the mean. Nature 189: 732–35.
#'
#' Hughes G, Madden LV. 1992. Aggregation and incidence of disease. Plant
#' Pathology 41 (6): 657–660.
#' \doi{10.1111/j.1365-3059.1992.tb02549.x}
#'
#' Madden LV, Hughes G, van den Bosch F. 2007. Spatial aspects of epidemics -
#' III: Patterns of plant disease. In: The study of plant disease epidemics,
#' 235–278. American Phytopathological Society, St Paul, MN.
#'
#' @export
#------------------------------------------------------------------------------#
power_law <- function(data, log_base = exp(1), ...) {

    call <- match.call()
    # Checks:
    stopifnot(is.list(data))
    if (length(data) < 2) {
        stop("Less than 2 points is not enough to perform linear regressions.")
    }
    object_class <- unique(vapply(data, function(x) class(x)[[1L]],
                                  character(1L)))
    stopifnot(length(object_class) == 1)
    stopifnot(object_class %in% c("count", "incidence"))

    # Perform power law analysis:
    switch(object_class,
           "count" = {
               name <- "Taylor's Power Law"
               data <- get_fmt_obs(data, type = object_class)
               x    <- vapply(data, mean, numeric(1L))
               y    <- vapply(data, var, numeric(1L))
           },
           "incidence" = {
               name <- "Binary Power Law"
               data <- get_fmt_obs(data, type = object_class)
               # For incidence data as proportions:
               # v_t = p(1 - p)/n (Madden & Hughes, 1995)
               x    <- vapply(data, function(obj) {
                   with(obj, (mean(p) * (1 - mean(p))) / mean(n))
               }, numeric(1L))
               y    <- vapply(data, function(obj) var(obj$p), numeric(1L))
           }
    )
    coord_obs <- data.frame(x = x, y = y)
    if (log_base == exp(1)) {
        # To get a nice display with print and summary
        model_formula <- as.formula(log(y) ~ log(x))
    } else {
        model_formula <- as.formula(bquote(
            log(y, base = .(log_base)) ~ log(x, base = .(log_base))
        ))
    }
    model     <- lm(model_formula, ...)
    y_the     <- predict(model, type = "response")
    coord_the <- data.frame(x = x, y = log_base^(y_the))

    # Retrieve summary matrice of coefficients, and eventually add some extra
    # estimates:
    param <- coef(summary(model))
    rownames(param) <- paste0(rownames(param), ": ", c("log_base(Ar)", "b"))
    switch (object_class,
        "count" = {
            # Nothing to do.
        },
        "incidence" = {
            n  <- mean(vapply(data, function(obj) mean(obj$n), numeric(1L))) ### PAS TOP
            new_param <- list(
                Ai = estimate_param(model, bquote(.(log_base)^x1)),
                ai = estimate_param(model, bquote(.(log_base)^x1 * .(n)^(-x2))),
                AI = estimate_param(model, bquote(.(log_base)^x1 * .(n)^(2 * (1 - x2)))),
                aI = estimate_param(model, bquote(.(log_base)^x1 * .(n)^(2 - x2)))
            )
            param <- rbind(param, do.call(rbind, lapply(new_param, unlist))) # TODO: Not clean
        }
    )

    # Return the following object:
    structure(list(call       = call, # TODO: Add more information about the transformation?
                   data_class = object_class,
                   name       = name,
                   data       = data, # TODO: Useful ??? (not in spatial_hier)
                   model      = model,
                   param      = param, # TODO: Where in n????
                   log_base   = log_base,
                   coord_obs  = coord_obs,
                   coord_the  = coord_the),
              class = "power_law")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.power_law <- function(x, ...) {
    cat(x$name, ":\n",
        "Power law analysis for '", x$data_class[1L], "' data.\n", sep = "")
    cat("\nCoefficients:\n")
    print(coef(x$model))
    cat("\n")
    invisible(x)
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
summary.power_law <- function(object, ...) {
    structure(list(call         = object$call,
                   model        = object$model,
                   coefficients = object$param
    ), class = "summary.power_law")
}

#------------------------------------------------------------------------------#
#' @method print summary.power_law
#' @export
#------------------------------------------------------------------------------#
print.summary.power_law <- function(x, ...) {
    res <- summary.lm(x$model, ...)
    res$call <- x$call
    res$coefficients <- x$coefficients
    print(res)
    invisible(res)
}

#------------------------------------------------------------------------------#
# Plot results of a power law analysis
#
# Plot results of a power law analysis.
#
# @param x A \code{\link{power_law}} object.
# @param scale Logarithmic or standard linear scale to display the results?
# @param observed Logical.
# @param model Logical.
# @param random Logical. Theoretical Random distribution.
# @param ... Additional arguments to be passed to other methods.
# TODO: Dashed lines indicate the cases where both variances are equal, which suggests an absence of aggregation.
# Points should lie on this line if ...
#
#' @export
#------------------------------------------------------------------------------#
plot.power_law <- function(x, ..., scale = c("logarithmic", "linear"),
                           observed = TRUE, model = TRUE, random = TRUE) {
    scale <- match.arg(scale)
    log_base <- x$log_base
    gg <- ggplot()

    switch (scale,
        "logarithmic" = {
            log_base_name <- ifelse(log_base == exp(1), "e",
                                    as.character(log_base))
            gg <- gg +
                labs(x = bquote(log[.(log_base_name)] * "(binomial variance)"),
                     y = bquote(log[.(log_base_name)] * "(observed variance)"))
            if (observed) {
                gg <- gg + geom_point(data = log(x$coord_obs, base = log_base),
                                      aes(x, y), ...)
            }
            if (model) {
                gg <- gg + geom_line(data  = log(x$coord_the, base = log_base),
                                     aes(x, y), ...)
            }
            if (random) {
                gg <- gg + geom_line(data  = log(x$coord_the, base = log_base),
                                     aes(x, x), linetype = "dashed", ...)
            }
        },
        "linear" = {
            gg <- gg + labs(x = "Binomial variance", y = "Observed variance")
            if (observed) {
                gg <- gg + geom_point(data = x$coord_obs, aes(x, y), ...)
            }
            if (model) {
                gg <- gg + geom_line(data  = x$coord_the, aes(x, y), ...)
            }
            if (random) {
                gg <- gg + geom_line(data  = x$coord_the, aes(x, x),
                                     linetype = "dashed", ...)
            }
        }
    )
    gg <- gg + theme_bw()
    print(gg)
    invisible(NULL)
}


#==============================================================================#
# a2a
#==============================================================================#

#------------------------------------------------------------------------------#
#' Easily switch between different power law formulations.
#'
#' \code{a2a} was designed to avoid headaches that are likely to occur when
#' working with different formulations of the binomial power law analysis.
#'
#' The binomial power law can be expressed as: \eqn{s_y^2 = (intercept)(s_{bin}^2)^b}.
#' But different forms of (intercept) are possible depending on the formulation of the
#' binomial power law.
#' \tabular{ccccc}{
#'       \tab Ai         \tab ai      \tab AI         \tab aI      \cr
#'    Ai \tab 1          \tab n^b     \tab n^(2(b-1)) \tab n^(b-2) \cr
#'    ai \tab n^(-b)     \tab 1       \tab n^(b-2)    \tab n^(-2)  \cr
#'    AI \tab n^(2(1-b)) \tab n^(2-b) \tab 1          \tab n^(-b)  \cr
#'    aI \tab n^(2-b)    \tab n^2     \tab n^b        \tab 1       \cr
#' }
#'
#' @param x Intercept parameter to be converted or a named list with the
#'     parameter to be converted ("Ai", "ai", "AI" or "aI"), the slope
#'     ("slope"), and the number of individual per sampling unit ("n").
#' @param from Kind of the input intercept parameter ("Ai", "ai", "AI" or "aI").
#' @param to Desired kind for the ouput intercept parameter ("Ai", "ai", "AI" or
#'     "aI").
#' @param slope Slope parameter.
#' @param n Number of individuals per sampling unit.
#' @param ... Additional arguments to be passed to other methods.
#'
#' @returns A numeric vector.
#'
#' @examples
#' # Values from the power_law() example:
#' Ai    <- 38.6245
#' slope <- 1.9356
#' n     <- 9
#'
#' # Usual function call syntax:
#' a2a(Ai, slope, n, from = "Ai", to = "ai")
#'
#' # Other syntaxes:
#' inputs <- list(Ai = Ai, slope = slope, n = n)
#' a2a(inputs, "ai")
#' require(magrittr)
#' inputs %>% a2a("ai")
#'
#' @export
#------------------------------------------------------------------------------#
a2a <- function(x, ...) UseMethod("a2a")

#------------------------------------------------------------------------------#
#' @keywords internal
#------------------------------------------------------------------------------#
a2a_internal <- function(intercept, b, n, from, to) {
    dico <- expand.grid(from = c("Ai", "ai", "AI", "aI"),
                        to = c("Ai", "ai", "AI", "aI"),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    #-----------------------------------------------------------------#
    #             | col Ai     |col ai  | col AI     | col aI
    dico$coef <- c(1,           n^b,     n^(2*(b-1)), n^(b-2), # row Ai
                   n^(-b),      1,       n^(b-2),     n^(-2),  # row ai
                   n^(2*(1-b)), n^(2-b), 1,           n^(-b),  # row AI
                   n^(2-b),     n^2,     n^b,         1)       # row aI
    #-----------------------------------------------------------------#
    item <- dico[which(dico$from == from & dico$to == to), ]
    res  <- intercept * item[["coef"]]
    attr(res, "param") <- c(intercept = intercept,
                            coef      = item[["coef"]],
                            slope     = b,
                            n         = n)
    attr(res, "class") <- c("a2a", "numeric")
    res
}

#------------------------------------------------------------------------------#
#' @rdname a2a
#' @export
#------------------------------------------------------------------------------#
a2a.numeric <- function(x, slope, n, # Here, x = intercept
                        from = c("Ai", "ai", "AI", "aI"),
                        to   = c("Ai", "ai", "AI", "aI"), ...) {

    # Checks and variable allocation:
    from      <- match.arg(from)
    to        <- match.arg(to)
    intercept <- x
    b         <- slope

    # Perform the conversion:
    a2a_internal(intercept, b, n, from, to)
}

#------------------------------------------------------------------------------#
#' @rdname a2a
#' @export
#------------------------------------------------------------------------------#
a2a.list <- function(x, # Here, x: list with 3 elements: intercept, slope and n.
                     to   = c("Ai", "ai", "AI", "aI"), ...) {

    # Checks and variable allocation:
    stopifnot(length(x) == 3)
    stopifnot(all(c("slope", "n") %in% names(x)))
    cases <- c("Ai", "ai", "AI", "aI")
    from  <- cases[cases %in% names(x)]
    stopifnot(length(from) == 1)
    to    <- match.arg(to)
    intercept <- x[[from]]
    b     <- x[["slope"]]
    n     <- x[["n"]]

    # Perform the conversion:
    a2a_internal(intercept, b, n, from, to)
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.a2a <- function(x, ...) cat(x, "\n", sep = "")

