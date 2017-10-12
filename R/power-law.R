#------------------------------------------------------------------------------#
#' @include intensity-classes.R
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' Taylor's and binary power laws
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
#' @param x A list of \code{intensity} objects (\code{count} or
#'     \code{incidence} objects).
#' @param log_base Logarithm base to be used.
#' @param ... Not yet implemented.
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
#' # given time at a given location.
#' my_power_law <- power_law(my_data)
#' my_power_law
#' summary(my_power_law)
#' plot(my_power_law)
#'
#' @references
#'
#' Taylor LR. 1961. Aggregation, variance and the mean. Nature 189: 732–35.
#'
#' Hughes G, Madden LV. 1992. Aggregation and incidence of disease. Plant
#' Pathology 41 (6): 657–660.
#' \href{http://dx.doi.org/10.1111/j.1365-3059.1992.tb02549.x}{doi:10.1111/j.1365-3059.1992.tb02549.x}
#'
#' Madden LV, Hughes G, van den Bosch F. 2007. Spatial aspects of epidemics -
#' III: Patterns of plant disease. In: The study of plant disease epidemics,
#' 235–78. American Phytopathological Society, St Paul, MN.
#'
#' @export
#------------------------------------------------------------------------------#
power_law <- function(list, log_base = exp(1), ...) {

    # Checks:
    stopifnot(is.list(list))
    if (length(list) < 2) {
        stop("Less than 2 points is not enough to perform linear regressions.")
    }
    object_class <- unique(vapply(list, function(x) class(x)[[1L]],
                                  character(1L)))
    stopifnot(length(object_class) == 1)
    stopifnot(object_class %in% c("count", "incidence"))

    # Perform power law analysis:
    switch(object_class,
           "count" = {
               data <- lapply(list, function(obj) {
                   data_all  <- map_data(obj)$r
                   data_noNA <- data_all[complete.cases(data_all)]
                   if (length(data_noNA) < length(data_all)) {
                       warning("Missing cases were dropped.")
                   }
                   data_noNA
               })
               x    <- vapply(data, function(obj) mean(obj), numeric(1L))
               y    <- vapply(data, function(obj) var(obj), numeric(1L))
           },
           "incidence" = {
               data <- lapply(list, function(obj) {
                   mapped_data <- map_data(obj)
                   data_all  <- data.frame(p = (mapped_data$r / mapped_data$n),
                                           n = mapped_data$n)
                   data_noNA <- data_all[complete.cases(data_all), ]
                   if (nrow(data_noNA) < nrow(data_all)) {
                       warning("Missing cases were dropped.")
                   }
                   data_noNA
               })
               # For incidence data as proportions:
               # v_t = p(1 - p)/n (Madden & Hughes, 1995)
               x    <- vapply(data, function(obj) {
                   with(obj, (mean(p) * (1 - mean(p))) / mean(n))
               }, numeric(1L))
               y    <- vapply(data, function(obj) var(obj$p), numeric(1L))
           }
    )
    coord_obs <- data.frame(x = x, y = y)
    model_formula <- as.formula(bquote(
        log(y, base = .(log_base)) ~ log(x, base = .(log_base))
    ))
    model     <- lm(model_formula, ...)
    y_the     <- predict(model, type = "response")
    coord_the <- data.frame(x = x, y = log_base^(y_the))

    # Retrieve summary matrice of coefficients, and eventually add some extra
    # estimates:
    par <- coef(summary(model))
    switch (object_class,
        "count" = {
            # Nothing to do.
        },
        "incidence" = {
            n  <- mean(vapply(data, function(obj) mean(obj$n), numeric(1L))) ### PAS TOP
            Ar <- estimateCoef(model, bquote(.(log_base)^x1))
            ar <- estimateCoef(model, bquote(.(log_base)^x1 * .(n)^(-x2)))
            AR <- estimateCoef(model, bquote(.(log_base)^x1 * .(n)^(2 * (1 - x2))))
            aR <- estimateCoef(model, bquote(.(log_base)^x1 * .(n)^(2 - x2)))
            par <- rbind(par, unlist(Ar), unlist(ar), unlist(AR), unlist(aR))
            rownames(par) <- c("log_base(Ar)", "b", "Ar", "ar", "AR", "aR")
        }
    )

    # Return the following object:
    structure(list(call      = match.call(),
                   data      = data,
                   model     = model,
                   par       = par,
                   log_base  = log_base,
                   coord_obs = coord_obs,
                   coord_the = coord_the),
              class = "power_law")
}

#------------------------------------------------------------------------------#
#' Plot results of a power law analysis
#'
#' Plot results of a power law analysis.
#'
#' @param x A \code{\link{power_law}} object.
#' @param scale Logarithmic or standard linear scale to display the results?
#' @param observed Logical.
#' @param model Logical.
#' @param random Logical. Theoretical Random distribution.
#' Dashed lines indicate the cases where both variances are equal, which suggests an absence of aggregation.
#' Points should lie on this line if ...
#'
#' @examples
#'
#' plot(my_power_law, scale = "log")
#' plot(my_power_law, scale = "lin")
#'
#' @export
#------------------------------------------------------------------------------#
plot.power_law <- function(x, ..., scale = c("logarithmic", "linear"), observed = TRUE,
                           model = TRUE, random = TRUE) {
    scale <- match.arg(scale)
    data_obs <- x$coord_obs
    data_the <- x$coord_the
    log_base <- x$log_base

    switch (scale,
        "logarithmic" = {
            log_base_name <- ifelse(log_base == exp(1), "e", as.character(log_base))
            gg <- list(
                labs(x = bquote(log[.(log_base_name)] * "(binomial variance)"),
                     y = bquote(log[.(log_base_name)] * "(observed variance)")),
                # Below, switch is a quick & convenient way to say that if
                # requirement is TRUE (i.e. = 1), then return the following
                # instruction. Otherwise if requirement is FALSE (i.e. = 0),
                # return NULL.
                switch(observed, geom_point(data = log(data_obs, base = log_base),
                                            aes(x, y), ...)),
                switch(model,    geom_line(data  = log(data_the, base = log_base),
                                           aes(x, y), ...)),
                switch(random,   geom_line(data  = log(data_the, base = log_base),
                                           aes(x, x), linetype = "dashed", ...)),
                theme_bw()
            )
        },
        "linear" = {
            gg <- list(
                labs(x = "Binomial variance", y = "Observed variance"),
                switch(observed, geom_point(data = data_obs, aes(x, y), ...)),
                switch(model,    geom_line(data  = data_the, aes(x, y), ...)),
                switch(random,   geom_line(data  = data_the, aes(x, x),
                                           linetype = "dashed", ...)),
                theme_bw()
            )
        }
    )
    ggplot() + gg
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.power_law <- function(x, ...) {
    cat("# Power law analysis:\n")
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(coef(x$model))
    cat("\n")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
summary.power_law <- function(object, ...) {
    # TODO: Bien regarder la structure d'un summary.lm pour bien comprendre
    # tous ses éléments.
    summary_model <- summary(object$model)
    summary_model$call <- object$call
    summary_model$coefficients <- object$par
    structure(summary_model, class = "summary.power_law")
}

#------------------------------------------------------------------------------#
#' @method print summary.power_law
#' @export
#------------------------------------------------------------------------------#
print.summary.power_law <- function(x, ...) stats:::print.summary.lm(x, ...)


#==============================================================================#
# a2a
#==============================================================================#

#------------------------------------------------------------------------------#
#' a2a: A wAy to pAinlessly switch between different power LAw formulAtions
#'
#' \code{a2a} was designed to avoid headaches that are likely to occur when
#' working with different formulations of the binomial power law analysis.
#'
#' The binomial power law can be expressed as: \eqn{s_y^2 = (intercept)(s_{bin}^2)^b}.
#' But different forms of (intercept) are possible depending on the formulation of the
#' binomial power law.
#' \tabular{ccccc}{
#'       \tab Ar         \tab ar      \tab AR         \tab aR      \cr
#'    Ar \tab 1          \tab n^b     \tab n^{2(b-1)} \tab n^{b-2} \cr
#'    ar \tab n^{-b}     \tab 1       \tab n^{b-2}    \tab n^{-2}  \cr
#'    AR \tab n^{2(1-b)} \tab n^{2-b} \tab 1          \tab n^{-b}  \cr
#'    aR \tab n^{2-b}    \tab n^2     \tab n^b        \tab 1       \cr
#' }
#'
#' @param x Intercept parameter to be converted or a named list with the
#'     parameter to be converted ("Ar", "ar", "AR" or "aR"), the slope
#'     ("slope"), and the number of individual per sampling unit ("n").
#' @param from Kind of the input intercept parameter ("Ar", "ar", "AR" or "aR").
#' @param to Desired kind for the ouput intercept parameter ("Ar", "ar", "AR" or
#'     "aR").
#' @param slope Slope parameter.
#' @param n Number of individuals per sampling unit.
#'
#' @examples
#' # Values from the power_law() example:
#' Ar    <- 38.6245
#' slope <- 1.9356
#' n     <- 9
#'
#' # Usual function call syntax:
#' a2a(Ar, slope, n, from = "Ar", to = "ar")
#'
#' # Other syntaxes:
#' inputs <- list(Ar = Ar, slope = slope, n = n)
#' a2a(inputs, "ar")
#' require(magrittr)
#' inputs %>% a2a("ar")
#'
#' @export
#------------------------------------------------------------------------------#
a2a <- function(x, ...) UseMethod("a2a")

#------------------------------------------------------------------------------#
#' @keywords internal
#------------------------------------------------------------------------------#
a2a_internal <- function(intercept, b, n, from, to) {
    dico <- expand.grid(from = c("Ar", "ar", "AR", "aR"),
                        to = c("Ar", "ar", "AR", "aR"),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    #-----------------------------------------------------------------#
    #             | col Ar     |col ar  | col AR     | col aR
    dico$coef <- c(1,           n^b,     n^(2*(b-1)), n^(b-2), # row Ar
                   n^(-b),      1,       n^(b-2),     n^(-2),  # row ar
                   n^(2*(1-b)), n^(2-b), 1,           n^(-b),  # row AR
                   n^(2-b),     n^2,     n^b,         1)       # row aR
    #-----------------------------------------------------------------#
    item <- dico[which(dico$from == from & dico$to == to), ]
    res  <- intercept * item[["coef"]]
    attr(res, "par") <- c(intercept = intercept,
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
                        from = c("Ar", "ar", "AR", "aR"),
                        to   = c("Ar", "ar", "AR", "aR")) {

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
                     to   = c("Ar", "ar", "AR", "aR")) {

    # Checks and variable allocation:
    stopifnot(length(x) == 3)
    stopifnot(all(c("slope", "n") %in% names(x)))
    cases <- c("Ar", "ar", "AR", "aR")
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

