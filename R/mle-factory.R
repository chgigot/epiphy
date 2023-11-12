#------------------------------------------------------------------------------#
#' @include utils.R
#' @include betabinom.R
#------------------------------------------------------------------------------#
NULL

#==============================================================================#
# Class smle
#==============================================================================#

#------------------------------------------------------------------------------#
#' Simple maximum likelihood estimation
#'
#' By default, this function performs a maximum likelihood estimation for one or
#' several parameters, but it can be used for any other optimization problems.
#' The interface is intented to be rather simple while allowing more advanced
#' parametrizations.
#'
#' The \code{\link[stats]{optim}} tool does the hard work under the hood. Extra
#' arguments (e.g. method of optimization to be used) can be passed to
#' \code{\link[stats]{optim}} through the \code{...} argument. Note that
#' contrary to the default \code{\link[stats]{optim}} arguments, \code{smle}
#' tries to solve a maximization problem using the method "L-BFGS-B" by default
#' (see \code{\link[stats]{optim}} documentation for more information).
#'
#' @param data The data set to work with. It can be a vector (if there is only
#'     one variable), a data frame (if there is one or more variables) or an
#'     \code{\link{intensity}} object.
#' @param f A function to be maximized, typically a log-likelihood function.
#'     This function must have only two arguments: \code{data} and \code{param},
#'     which must correspond to the \code{data} argument of \code{smle} and a
#'     named vector of the parameter(s) to be estimated.
#' @param param_init Either a named vector with proposed initial values of the
#'     parameter(s), or a function that returns such a vector. This parameter
#'     is not needed if the parameter \code{param} of \code{f} is already
#'     provided with such a named vector.
#' @param max Does \code{f} need to be maximized? Set to \code{FALSE} to require
#'     a minimization of \code{f}.
#' @param ... Additional arguments to be passed to \code{\link[stats]{optim}}.
#'
#' @returns
#' An object of class \code{smle}, which is a list containing the following
#' components:
#' \tabular{ll}{
#'     \code{call}                    \tab The call. \cr
#'     \code{coef}                    \tab The estimated coefficients. \cr
#'     \code{coef_se}                 \tab The standard errors of the estimated coefficients. \cr
#'     \code{vcov}                    \tab The variance-covariance matrix of the estimated coefficients. \cr
#'     \code{data}                    \tab The \code{data} parameter. \cr
#'     \code{f}                       \tab The \code{f} parameter. \cr
#'     \code{nobs}                    \tab The number of observations. \cr
#'     \code{full_input, full_output} \tab The full input and output of the \code{optim} function. \cr
#' }
#'
#' @examples
#' set.seed(12345)
#' data <- rlogis(100, location = 5, scale = 2)
#' ll_logis <- function(data, param = c(location = 0, scale = 1)) {
#'     sum(dlogis(x = data, location = param[["location"]],
#'                scale = param[["scale"]], log = TRUE))
#' }
#' res <- smle(data, ll_logis)
#' res
#' summary(res)
#'
#' # Using the magrittr syntax:
#' require(magrittr)
#' data %>% smle(ll_logis)
#'
#' # Comparision with the output of fitdistr (MASS package), which works for a
#' # limited number of predefined distributions:
#' require(MASS)
#' fitdistr(data, "logistic")
#'
#' # Example with an intensity object:
#' require(magrittr)
#' require(dplyr)
#' data <- tomato_tswv$field_1929 %>%
#'     filter(t == 1) %>%
#'     incidence() %>%
#'     clump(unit_size = c(x = 3, y = 3))
#' ll_betabinom <- function(data, param) {
#'     sum(dbetabinom(x = data[["i"]], size = data[["n"]],
#'                    prob = param[["prob"]], theta = param[["theta"]],
#'                    log = TRUE))
#' }
#' epsilon <- 1e-7
#' res <- smle(data, ll_betabinom, param_init = c(prob = 0.5, theta = 0.5),
#'             lower = c(prob  = 0 + epsilon,
#'                       theta = 0 + epsilon),
#'             upper = c(prob = 1 - epsilon,
#'                       theta = Inf))
#' res
#' summary(res)
#'
#' param_init <- data.frame(lower = c(0 + epsilon, 0 + epsilon),
#'                          start = c(0.5, 0.5),
#'                          upper = c(1 - epsilon, Inf))
#' rownames(param_init) <- c("prob", "theta")
#' res <- smle(data, ll_betabinom, param_init)
#' res
#' summary(res)
#'
#' @keywords internal
#' @export
#------------------------------------------------------------------------------#
smle <- function(data, ...) UseMethod("smle")

#------------------------------------------------------------------------------#
#' @rdname smle
#' @export
#------------------------------------------------------------------------------#
smle.default <- function(data, f, param_init, max = TRUE, ...) {
    call <- match.call() # Vérifier que c'est OK lorsque l'on travaille avec un obj intensity
    # Checks and data preparation: TODO
    dots <- list(...)
    if (!is.vector(data) && !is.data.frame(data)) {
        stop("data must be a vector or a data frame.")
    }
    data <- na.omit(data)
    fullcoef <- formals(f)
    if (!missing(param_init)) {
        if (is.numeric(param_init)) {
            # STOP if not a named vector
        } else {
            # TODO: Below what if it's not a data frame or a matrix (if we have just start)?
            if (is.function(param_init)) {
                param_name <- get_param_name_from_body(f)
                param_init <- param_init(data = data, name = param_name) # function param_init is overwritten here
                # Checks... TODO
                ## browser()
            }
            if (!is.data.frame(param_init)) {
                param_init <- as.data.frame(param_init)
            }
            get_named_col <- function(x, name) {
                if (!is.null(col <- x[[name]])) setNames(x[[name]], rownames(x))
                else NULL # Works, return NULL in a list => nothing happen
            }
            if (is.null(dots$lower)) dots$lower <- get_named_col(param_init, "lower")
            if (is.null(dots$upper)) dots$upper <- get_named_col(param_init, "upper")
            param_init <- get_named_col(param_init, "start") # data frame param_init is overwritten with a vector here
        }
    } else {
        # If param is a call (named vector, function):
        if (is.call(fullcoef$param)) {
            param_init <- eval(fullcoef$param)
        } else {
            stop("No provided initial values for the parameter(s).")
        }
    }
    # To perform a maximization problem (instead of the default minimization
    # problem in optim), fnscale must be < 0:
    #if (max) {
    #    control <- list(fnscale = -1)
    #} else {
    #    control <- list()
    #}
    #if (!is.null(dots$control)) {
    #    control <- modifyList(control, dots$control)
    #}
    if (is.null(dots$control$fnscale)) {
        if (max) {
            dots$control$fnscale <- -1
        }
    }
    if (is.null(dots$method)) {
        dots$method <- "L-BFGS-B"
    }
    if (is.null(dots$lower)) dots$lower <- eval(formals(optim)[["lower"]])
    if (is.null(dots$upper)) dots$upper <- eval(formals(optim)[["upper"]])

    # Add missing optim's arguments
    dots <- c(dots, list(
        fn      = f,
        data    = data,
        par     = param_init,
        hessian = TRUE
    ))
    # Perform the optimization itself:
    out <- NULL
    try(out <- do.call(optim, dots), silent = TRUE)
    if (!is.null(out)) {
        coef    <- out$par
        # vcov: Estimator of the asymptotic covariance matrix (Fisher info,
        # "observed information"). We need to use solve(-hessian) with a
        # maximization problem, and solve(hessian) with a minimization problem.
        # ref: https://stats.stackexchange.com/questions/27033/in-r-given-an-output-from-optim-with-a-hessian-matrix-how-to-calculate-paramet
        if (dots$control$fnscale < 0) {
            vcov <- solve(-out$hessian)
        } else {
            vcov <- solve(out$hessian)
        }
        coef_se <- sqrt(diag(vcov))
    } else {
        ## TODO: BEG TMP
        ##dots$method  <- "Nelder-Mead"
        ##dots$hessian <- FALSE
        ##try(out <- do.call(optim, dots), silent = TRUE)
        ## ref: https://stackoverflow.com/questions/28185387/non-finite-finite-difference-value-many-data-become-inf-and-na-after-exponentia
        ## TODO: END TMP
        coef    <- NULL
        vcov    <- NULL
        coef_se <- NULL
    }

    structure(list(call = call, coef = coef, coef_se = coef_se, vcov = vcov,
                   data = data, f = f, nobs = NROW(data),
                   full_input = dots, full_output = out),
              class = "smle")
}

#------------------------------------------------------------------------------#
#' @rdname smle
#' @export
#------------------------------------------------------------------------------#
smle.intensity <- function(data, f, param_init, max = TRUE, ...) {
    call <- match.call()
    mapped_data <- map_data(data)
    res <- smle.default(data = mapped_data, f = f, param_init = param_init,
                        max = max, ...)
    res$call <- call
    res
}

#------------------------------------------------------------------------------#
#' @inherit stats::coef title description return
#' @inheritParams stats::coef
#' @keywords internal
#' @export
#------------------------------------------------------------------------------#
coef.smle <- function(object, ...) object$coef

#------------------------------------------------------------------------------#
#' @inherit stats::vcov title description return
#' @inheritParams stats::vcov
#' @keywords internal
#' @export
#------------------------------------------------------------------------------#
vcov.smle <- function(object, ...) object$vcov

#------------------------------------------------------------------------------#
#' Extract log-likelihood
#'
#' This function returns the maximal log-likelihood estimated with
#' \code{\link{smle}}, if \code{f} returned log-likelihood value.
#'
#' @inherit stats::logLik return
#' @inheritParams stats::logLik
#' @keywords internal
#' @export
#------------------------------------------------------------------------------#
logLik.smle <- function(object, ...) {
    # To get the (log-)likelihood, do not forget the minus when it was
    # a minimization problem (optim default):
    val <- -object$full_output$value
    if (!is.null(fnscale <- object$full_input$control$fnscale)) {
        if (fnscale < 0) {
            val <- object$full_output$value
        }
    }
    attr(val, "df")   <- length(object$coef) # To double check (bbmle et stats4)
    attr(val, "nobs") <- object$nobs
    structure(val, class = "logLik")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.smle <- function(x, ...) {
    if (is.null(x$full_output)) {
        cat("Maximum likelihood procedure did not succeed.\n")
    } else {
        invisible(lapply(seq_len(length(x$coef)), function(i1) {
            cat(names(coef(x))[i1], ": ", formatC(coef(x)[i1]),
                " (\u00b1 ", formatC(x$coef_se[i1]), ")\n", sep = "")
        }))
    }
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
summary.smle <- function(object, ...) {
    cmat <- cbind(coef(object), object$coef_se)
    cmat <- cbind(cmat, cmat[, 1] / cmat[, 2])
    cmat <- cbind(cmat, 2 * pnorm(abs(cmat[, 3]), lower.tail = FALSE))
    colnames(cmat) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
    res <- list(call_orig = object$call,
                cmat      = cmat)
    structure(res, class = "summary.smle")
}

#------------------------------------------------------------------------------#
#' @method print summary.smle
#' @export
#------------------------------------------------------------------------------#
print.summary.smle <- function(x, ...) {
    cat("Maximum likelihood estimation\n\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    printCoefmat(x$cmat)
    # cat("\n-2 log L:", x@m2logL, "\n")
}

#------------------------------------------------------------------------------#
#' @method coef summary.smle
#' @export
#------------------------------------------------------------------------------#
coef.summary.smle <- function(object, ...) object$cmat



#==============================================================================#
# Method of Moments Estimation
# Used to find initial estimates of the parameters.
#==============================================================================#
mme_pois <- function(data, name, bounds = TRUE) {
    data <- get_std_named_df(data, "i")
    x <- data[["i"]]
    m1 <- mean(x) # na.rm = TRUE?
    fmt_init(data, name, bounds = bounds,
             lambda = quote(c(0 + epiphy_env$epsilon,
                              m1,
                              Inf)))
}

mme_nbinom <- function(data, name, bounds = TRUE) {
    data <- get_std_named_df(data, "i")
    x <- data[["i"]]
    m1 <- mean(x)   # na.rm = TRUE?
    m2 <- mean(x^2) # na.rm = TRUE?
    s2 <- var(x)    # na.rm = TRUE?
    fmt_init(data, name, bounds = bounds,
             mu = quote(c(0 + epiphy_env$epsilon,
                          m1,
                          Inf)),
             k  = quote(c(0 + epiphy_env$epsilon,
                          m1^2/(s2 - m1), # Madden et al, 2007, p. 249
                          Inf)))
}

mme_binom <- function(data, name, bounds = TRUE) {
    data <- get_std_named_df(data, c("i", "n"))
    x <- data[["i"]]
    n <- data[["n"]][[1]] # if pas de n, alors n = max(x) /// et c'est pas tres propre
    m1 <- mean(x) # na.rm = TRUE?
    fmt_init(data, name, bounds = bounds,
             prob = quote(c(0 + epiphy_env$epsilon,
                            m1/n,
                            1 - epiphy_env$epsilon)))
}

mme_betabinom <- function(data, name, bounds = TRUE) { # TODO: To double-check
    data <- get_std_named_df(data, c("i", "n"))
    x <- data[["i"]]
    n <- data[["n"]][[1]] # if pas de n, alors n = max(x) /// et c'est pas tres propre
    m1 <- mean(x) # na.rm = TRUE?
    m2 <- mean(x^2) # na.rm = TRUE?
    s2 <- var(x) # x != proportion, but number of diseased individual per sampling unit.
                 # na.rm = TRUE?
    fmt_init(data, name, bounds = bounds,
             prob  = quote(c(0 + epiphy_env$epsilon,
                             m1/n,
                             1 - epiphy_env$epsilon)),
             theta = quote(c(0 + epiphy_env$epsilon,
                             (s2 - m1 * (1 - m1/n))/(n^2 * m1/n * (1 - m1/n) - s2),
                             Inf)),
             alpha = quote(c(0 + epiphy_env$epsilon,
                             (n * m1 - m2)/(n * (m2/m1 - m1 - 1) + m1),
                             Inf)),
             beta  = quote(c(0 + epiphy_env$epsilon,
                             (n - m1) * (n - m2/m1) / (n * (m2/m1 - m1 - 1) + m1),
                             Inf)))
}


#==============================================================================#
# Log-Likelihood functions
#==============================================================================#
ll_pois <- function(data, param) {
    data <- get_std_named_df(data, "i")
    sum(dpois(x = data[["i"]], lambda = param["lambda"], log = TRUE))
}

ll_nbinom <- function(data, param) {
    data <- get_std_named_df(data, "i")
    sum(dnbinom(x = data[["i"]], mu = param["mu"], size = param["k"],
                log = TRUE))
}

ll_binom <- function(data, param) {
    data <- get_std_named_df(data, c("i", "n"))
    sum(dbinom(x = data[["i"]], size = data[["n"]], prob = param[["prob"]],
               log = TRUE))
}

ll_betabinom <- function(data, param) {
    data <- get_std_named_df(data, c("i", "n"))
    sum(dbetabinom(x = data[["i"]], size = data[["n"]],
                   prob = param[["prob"]], theta = param[["theta"]],
                   #shape1 = param[["alpha"]], shape2 = param[["beta"]],
                   log = TRUE))
}


#==============================================================================#
# Maximum Likelihood Estimation
# Functions specific to a distribution
#==============================================================================#

#------------------------------------------------------------------------------#
#' Wrappers using maximum likelihood estimation for some distributions
#'
#' These functions are the core of the fitting processes performed in
#' \code{\link{fit_two_distr}}.
#'
#' @param data The data set to work with. It can be a vector (if there is only
#'     one variable), a data frame (if there is one or more variables) or an
#'     \code{\link{intensity}} object.
#'
#' @returns See \code{\link{smle}}
#'
#' @examples
#' set.seed(12345)
#' data <- rpois(100, lambda = 5)
#' res <- smle_pois(data)
#' res
#' summary(res)
#'
#' data <- count(aphids)
#' res <- smle_pois(data)
#' res
#' summary(res)
#'
#' @keywords internal
#' @name smle_wrappers
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' @rdname smle_wrappers
#' @export
#------------------------------------------------------------------------------
smle_pois <- structure(function(data) smle(data, ll_pois, mme_pois),
                       name = "Poisson")
# TODO: That would be great if we could use smle_*(data) instead of
# sme_*(data.frame(i = data)) if it's only a vector.

#------------------------------------------------------------------------------#
#' @rdname smle_wrappers
#' @export
#------------------------------------------------------------------------------
smle_nbinom <- structure(function(data) smle(data, ll_nbinom, mme_nbinom),
                         name = "Negative binomial")

#------------------------------------------------------------------------------#
#' @rdname smle_wrappers
#' @export
#------------------------------------------------------------------------------
smle_binom <- structure(function(data) smle(data, ll_binom, mme_binom),
                        name = "Binomial")

#------------------------------------------------------------------------------#
#' @rdname smle_wrappers
#' @export
#------------------------------------------------------------------------------#
smle_betabinom <- structure(function(data) smle(data, ll_betabinom, mme_betabinom),
                            name = "Beta-binomial")




