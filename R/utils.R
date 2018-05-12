#------------------------------------------------------------------------------#
# Utilities
#' @keywords internal
#------------------------------------------------------------------------------#
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    return(abs(x - round(x)) < tol)
}

#------------------------------------------------------------------------------#
#' @keywords internal
#------------------------------------------------------------------------------#
epiphy_env <- new.env()
epiphy_env$epsilon <- 1e-7

## Below useful for clumped_data in clump (= a list)?drop
#' @export
droplevels.list <- function(x, except = NULL, ...) {
    ix <- vapply(x, is.factor, logical(1L))
    if (!is.null(except)) {
        ix[except] <- FALSE
    }
    x[ix] <- lapply(x[ix], droplevels)
    x
}

#------------------------------------------------------------------------------#
#' Some link functions.
#'
#' Logit, probit and cloglog functions are available.
#' The logit and the logistic (with rev = TRUE), i.e. the inverse-logit functions.
#' Probit is a wrapper around \code{qnorm} (for \eqn{probit}) and \code{pnorm} (for \eqn{probit^{-1}})
#' Complementary log-log transformation.
#'
#' @param x A numeric vector.
#' @param rev The inverse of the function?
#'
#' @name link
#' @export
#------------------------------------------------------------------------------#
logit <- function(x, rev = FALSE) {
    if (!rev) log(x / (1 - x))
    else      1 / (1 + exp(-x))
}

#------------------------------------------------------------------------------#
#' @rdname link
#' @export
#------------------------------------------------------------------------------#
probit <- function(x, rev = FALSE) {
    if (!rev) qnorm(x)
    else      pnorm(x)
}

#------------------------------------------------------------------------------#
#' @rdname link
#' @export
#------------------------------------------------------------------------------#
cloglog <- function(x, rev = FALSE) {
    if (!rev) log(-log(1 - x))
    else      1 - exp(-exp(x))
}

#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Log-likelihood ratio test / likelihood ratio test statistic (LRS)
#
# idea from  lmtest::lrtest()
#
# Performs a Log-likelihood ratio test
#
# retravailler cette fonction pour qu'elle exploite directemenent lmtest::lrtest()
#
# random : random model
# aggregated : aggregated model
#
# not exported currently
#------------------------------------------------------------------------------#
llrtest <- function(random, aggregated) {
    llval_random     <- logLik(random)
    llval_aggregated <- logLik(aggregated)
    if (llval_aggregated < llval_random) {
        warning(paste0("logLik value is lower for aggregated model than for ",
                       "random model.\nCheck function arguments."))
    }
    value <- 2 * (llval_aggregated[[1]] - llval_random[[1]]) # [[1]] (or [1]) to only keep the ll value (not attributes)
    df    <- attr(llval_aggregated, "df") - attr(llval_random, "df") # Should be 1
    pval  <- pchisq(value, df, lower.tail = FALSE)
    chisq <- qchisq(pval, df)

    tab <- matrix(rep(NA, 8), nrow = 2)
    dimnames(tab) <- list(c("random :", "aggregated :"),
                          c("LogLik", "Df", "Chisq", "Pr(>Chisq)"))
    tab[, 1]  <- c(llval_random, llval_aggregated)
    tab[2, 2] <- df
    tab[2, 3] <- value
    tab[2, 4] <- pval
    title <- "Likelihood ratio test\n"
    structure(as.data.frame(tab),
              heading = title, # heading must be an attribute.
              class = c("anova", "data.frame"))
}


#------------------------------------------------------------------------------#
# Notes
#------------------------------------------------------------------------------#
# exp((AICmin−AICi)/2) peut être compris comme la probabilité pour que le ième
# candidat modèle minimise l'estimation de la perte d'information
# ref:Burnham et Anderson 2002, §6.4.5
#exp((118.9-152.1)/2)
#------------------------------------------------------------------------------#
# Christensen (1990 *) said: " The likelihood ratio test statistic is
# G2 = −2[log L(m0) − log L( m)], where ˆm0 is the MLE of m under the assumption
# that H0 is true and m is the MLE under the “unrestricted” model ".
# * Christensen, R. (1990)Log-Linear Models and Logistic Regression.
# Springer-Verlag New York, Inc., pp. 332-336.
#------------------------------------------------------------------------------#
# http://stackoverflow.com/questions/1826519/function-returning-more-than-one-value
# https://stat.ethz.ch/pipermail/r-help/2004-June/053343.html
# From: G. Grothendieck

# Some refs:
# - About negbinom (tutorial):
# http://rstudio-pubs-static.s3.amazonaws.com/2516_d204dc109a2c44d58899afb418ae3885.html
# http://data.princeton.edu/wws509/r/overdispersion.html
# http://stats.stackexchange.com/questions/10419/what-is-theta-in-a-negative-binomial-regression-fitted-with-r
# - About Hessian and Co.
# http://stats.stackexchange.com/questions/27033/in-r-given-an-output-from-optim-with-a-hessian-matrix-how-to-calculate-paramet
# http://stats.stackexchange.com/questions/68080/basic-question-about-fisher-information-matrix-and-relationship-to-hessian-and-s




#------------------------------------------------------------------------------#
# Estimate extra coef, based on deltamethod
# @export
#------------------------------------------------------------------------------#
estimate_param <- function(model, expr, type = c("t", "norm"), as_list = FALSE) { # as_list... NOT USED. List uniquement pour estimate mean / trouver autre nom
    # NEW !!!
    type <- match.arg(type)
    #expr_call <- substitute(expr) ## TODO: Non, parce que dans ce cas,on retrouve le call passé via expr (e.g. ".bquote(...)")
    formula_call <- deparse(expr, width.cutoff = 500)
    formula_call <- as.formula(paste0("~", formula_call))
    # END NEW !!!
    # resFormula <- as.formula(paste0("~", deparse(expr)))
    coefs <- coef(model)
    envir <- list2env(setNames(as.list(coefs), paste0("x", 1:length(coefs))))

    res      <- list()
    res$est  <- eval(expr, envir)
    res$se   <- msm::deltamethod(formula_call, coef(model), vcov(model))
    res$tval <- res$est / res$se
    res$pval <- ifelse((type == "t"), ## faire une condition plus intelligente (t and norm)
                         2 * pt(abs(res$tval), df.residual(model), lower.tail = FALSE), # TODO: df.residual with smle???
                         2 * pnorm(abs(res$tval), lower.tail = FALSE))

    # (1) Renvoyer les bons noms : par exemple:
    # Estimate Std. Error   z value        Pr(z)
    # ou
    # Estimate Std. Error t value Pr(>|t|)

    # (2) Renvoyer sous forme de matrix si demandé

    # (3) Vectorialisé si plusieurs paramètres demandés,
    # eg. estimate_param(model, list(Ar = bquote(.(log_base)^x1),
    #                                ar = bquote(.(log_base)^x1 * .(n)^(-x2))))

    res
}

estimate_param_ <- function(model, expr_name, type = c("t", "norm")) {
    type <- match.arg(type)
    expr <- as.name(expr_name)
    browser()
    estimate_param(model, expr, type)
}


# mat: existing mat of parameters
# new_param: list of new params
rbind_param <- function(base_mat, new_param) {
    rbind(base_mat, do.call(rbind, lapply(new_param, unlist)))
}


# @export
# extraCoef <- function(model) { ## TODO: USEFUL ???
#     #myfunc <- function(v1) {
#     #    deparse(substitute(v1))
#     #}
#     #
#     #myfunc(foo)
#     #[1] "foo"
#
#     # Retrieve result matrice (to which we will add extra estimates)
#     coefs   <- coef(stats4::summary(model)) # NON pas normal de devoir mettre bbmle ici, pas besoin de le faire dans estimate_param
#     # not need to do that forcoef(model) ... why?
#
#     # pour betabinom uniquement
#     extraCoefs = list(
#         p     = quote(x1 / (x1 + x2)),
#         theta = quote(1  / (x1 + x2)),
#         rho   = quote(1  / (x1 + x2 + 1))
#     )
#
#     newCoefs <- do.call(rbind, lapply(extraCoefs, function(i1) {
#         unlist(estimate_param(model, i1, type = "norm")) # faire plus intelligent ici avec le "norm"
#     }))
#     rbind(coefs, newCoefs)
# }

#------------------------------------------------------------------------------#
# Following function used in power_law and spatial-hier
# Get formatted observational variables
get_fmt_obs <- function(list, type = c("count", "incidence")) {
    type <- match.arg(type)
    switch (type,
        "count" = {
            lapply(list, function(obj) {
                data_all  <- map_data(obj)[["i"]]
                data_noNA <- data_all[complete.cases(data_all)]
                if (length(data_noNA) < length(data_all)) {
                    warning("Missing cases were dropped.")
                }
                data_noNA
            })
        },
        "incidence" = {
            lapply(list, function(obj) {
                mapped_data <- map_data(obj)
                data_all    <- data.frame(p = mapped_data[["i"]] / mapped_data[["n"]],
                                          n = mapped_data[["n"]])
                data_noNA   <- data_all[complete.cases(data_all), ]
                if (nrow(data_noNA) < nrow(data_all)) {
                    warning("Missing cases were dropped.")
                }
                data_noNA
            })
        }
    )
}


# utils used in mle-factory.R:
get_param_name_from_body <- function(f) {
    regex <- "param[\\[]{1,2}\"(.+?)\"[\\]]" # .+?: interogation mark: mean non-greedy behavior (Perl-compatible)
    f <- deparse(body(f))
    m <- gregexpr(regex, f, perl = TRUE)
    s <- regmatches(f, m)
    unique(unlist(lapply(s, function(ss) {
        if (length(ss) > 0) {
            gsub(regex, "\\1", ss, perl = TRUE)
        }
    })))
}

# always lower, start and upper,
# need documentation
fmt_init <- function(data, name, ..., bounds) {
    calling_env <- parent.frame()
    dots <- list(...)
    res <- do.call(rbind, lapply(setNames(name, name), function(val) {
        if (!is.null(call <- dots[[val]])) {
            names(call) <- c("", "lower", "start", "upper")
            res <- eval(call, envir = calling_env)
            # Checks:
            if (res["start"] < res["lower"]) res["start"] <- res["lower"]
            if (res["start"] > res["upper"]) res["start"] <- res["upper"]
            if (!bounds) res["start"] else res
        }
    }))
    as.data.frame(res)
}

tocamel <- function(x) {
    # source: https://stackoverflow.com/questions/6364783/capitalize-the-first-letter-of-both-words-in-a-two-word-string
    gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", x, perl=TRUE)
}


#------------------------------------------------------------------------------#
# TODO: to move to somewhere else?
#------------------------------------------------------------------------------#
chisq.test2 <- function(x, p, n_est, df, rescale.p = FALSE, ...) {

    if (missing(n_est) && missing(df)) {
        res <- stats::chisq.test(x = x, p = p, rescale.p = rescale.p, ...)
        res$data.name <- deparse(substitute(x))
        res
    } else {
        if (!missing(n_est) && !missing(df)) {
            stop("n_est or df must be provided, not both at the same time.")
        }
        if (!missing(n_est)) df <- length(x) - n_est - 1
        n <- sum(x)
        if (abs(sum(p) - 1) > sqrt(.Machine$double.eps)) {
            if (rescale.p) p <- p/sum(p)
            else stop("Probabilities must sum to 1.")
        }
        E <- n * p
        V <- n * p * (1 - p)
        statistic <- sum((x - E)^2/E)
        pVal <- pchisq(statistic, df = df, lower.tail = FALSE)
        names(statistic) <- "X-squared"
        names(df) <- "df"
        names(E) <- names(x)
        if (any(E < 5)) warning("Chi-squared approximation may be incorrect.")
        structure(list(statistic = statistic,
                       parameter = df,
                       p.value = pVal,
                       method = "Chi-squared goodness-of-fit with intrinsic null hypothesis",
                       data.name = deparse(substitute(x)),
                       observed = x,
                       expected = E,
                       residuals = (x - E)/sqrt(E),
                       stdres = (x - E)/sqrt(V)),
                  class = "htest")
    }
}




#------------------------------------------------------------------------------#
# Below : used in mle_factory
# data = vector, data.Frame or matrix
# order is important!!
get_std_named_df <- function(data, name) {
    if (is.data.frame(data) && all(name %in% colnames(data))) {
        # If it's already a std data frame (with potentially more columns that
        # only in name)
        data
    } else {
        stopifnot(all(sapply(data, is.numeric)))
        stopifnot(is.character(name))
        stopifnot(NCOL(data) == length(name))
        # Coerce to data.frame if it is a vector for instance
        data <- data.frame(data)
        setNames(data, name[1: ncol(data)])
    }
}


#==============================================================================#

#------------------------------------------------------------------------------#
#
#------------------------------------------------------------------------------#
as.long.data.frame <- function(x, ...) UseMethod("as.long.data.frame")


#------------------------------------------------------------------------------#
#' @method as.long.data.frame matrix
#------------------------------------------------------------------------------#
as.long.data.frame.matrix <- function(x, col_names = c("x", "y", "z"), ...) {
    id <- c(TRUE, TRUE)
    if (!is.null(dimnames(x))) {
        id <- vapply(dimnames(x), is.null, logical(1L))
    }
    dimnames(x)[id] <- list(seq_len(dim(x)[1L]), seq_len(dim(x)[2L]))[id]
    x <- as.data.frame(as.table(x))
    x[1:2] <- lapply(x[1:2], function(y) as.numeric(levels(y))[y])
    colnames(x) <- col_names
    x
}



#==============================================================================#


#------------------------------------------------------------------------------#
#' Retrieve vector or array indices
#'
#' \code{ind2sub} is a synonym for \code{\link[base]{arrayInd}} in \code{base} package.
#'
#' \code{ind2sub} is just an alias for \code{\link[base]{arrayInd}}
#'
#' @param ind Vector indices.
#' @param sub Array/matrix indices.
#'
#' @examples
#' set.seed(12345)
#' mat <- matrix(round(runif(6, min = 0, max = 10)), nrow = 2, ncol = 3)
#' ind2sub(4, dim(mat))
#' sub2ind(c(2, 2), dim(mat))
#' subs <- as.matrix(expand.grid(1:2,2:3))
#' sub2ind(subs, dim(mat))
#'
#' @keywords internal
#' @name indAndSub
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' @rdname indAndSub
#' @inheritParams base::arrayInd
#' @export
#------------------------------------------------------------------------------#
ind2sub <- base::arrayInd
#ind2sub <- function(ind, .dim, .dimnames = NULL, useNames = FALSE)
#    base::arrayInd(ind, .dim, .dimnames, useNames)
# http://stackoverflow.com/questions/4452039/rs-equivalent-to-ind2sub-sub2ind-in-matlab
# For ?-dimensional matrices
# Inspired by:
# http://tipstrickshowtos.blogspot.com/2011/09/fast-replacement-for-ind2sub.html
#ind2sub <- function(...) {
#    r = rem(idx-1,nrows)+1;
#    c = (idx-r)/nrows + 1;
#    # Equivalent to: [r,c] = ind2sub([nrows ncols],idx);
#}

#------------------------------------------------------------------------------#
#' @rdname indAndSub
#' @inheritParams base::arrayInd
#' @export
#------------------------------------------------------------------------------#
# For N-dimensional matrices
# Inspired by:
# http://tipstrickshowtos.blogspot.com/2010/02/fast-replacement-for-sub2ind.html
sub2ind <- function(sub, .dim) {
    if (!is.matrix(sub)) {
        sub <- matrix(sub, nrow = 1L)
    }
    stopifnot(ncol(sub) == length(.dim))
    nr <- nrow(sub)
    if (any(sub < 1L) || (any(sub > matrix(rep(.dim, each = nr), nrow = nr)))) {
        stop("Error: subscript out of bounds.")
    }
    ind <- apply(sub, 1L, function(coord) {
        sum(vapply(seq_len(length(coord)), function(k) {
            if (k == 1L) x <- coord[k]
            else         x <- (coord[k] - 1L) * prod(.dim[1L:(k-1L)])
            x
        }, numeric(1L)))
    })
    storage.mode(ind) <- "integer"
    ind
}












