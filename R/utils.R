#------------------------------------------------------------------------------#
# Utilities
#' @keywords internal
#------------------------------------------------------------------------------#
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    return(abs(x - round(x)) < tol)
}

#------------------------------------------------------------------------------#
#' Some link functions
#'
#' Logit, probit and cloglog functions are available.
#' The logit and the logistic (with rev = TRUE), i.e. the inverse-logit functions.
#' Probit is a wrapper around \code{qnorm} (for \eqn{probit}) and \code{pnorm} (for \eqn{probit^{-1}})
#' Complementary log-log transformation.
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
# not exported currently
#------------------------------------------------------------------------------#
llrtest <- function(modRand, modAgg) {
    llModRand <- as.numeric(logLik(modRand))
    llModAgg  <- as.numeric(logLik(modAgg))
    value <- 2 * (llModAgg - llModRand)
    df <- 1 #2 * (N - 1)
    pvalue <- 1 - pchisq(value, df) # use that : "lower.tail = FALSE)" instead
    chisq <- qchisq(pvalue, df)

    tab <- matrix(rep(NA, 8), nrow = 2)
    dimnames(tab) <- list(c("random :", "aggregated :"),
                          c("LogLik", "Df", "Chisq", "Pr(>Chisq)"))
    tab[, 1]  <- c(llModRand, llModAgg)
    tab[2, 2] <- df
    tab[2, 3] <- value
    tab[2, 4] <- pvalue

    title <- "Likelihood ratio test\n"
    #topnote <- paste0("Model (random)\nModel (agg)")
    #heading = c(title, topnote)
    structure(as.data.frame(tab), heading = title,
              class = c("anova", "data.frame"))
}


#------------------------------------------------------------------------------#
# Estimate extra coef, based on deltamethod
#' @export
#------------------------------------------------------------------------------#
estimateCoef <- function(model, expr, envir = NULL, type = c("t", "norm")) { # List uniquement pour estimate mean / trouver autre nom
    # NEW !!!
    type <- match.arg(type)
    resFormula <- paste0(deparse(expr), collapse = "") # Utilise si très longue formule !!!!!, car expr sur plusieurs lignes si longue formule... pas top
    resFormula <- as.formula(paste0("~", resFormula))
    # END NEW !!!
    # resFormula <- as.formula(paste0("~", deparse(expr)))
    res        <- list()
    if (is.null(envir)) {
        coefs <- coef(model)
        envir <- setNames(as.list(coefs), paste0("x", 1:length(coefs)))
    }
    res$est    <- eval(expr, envir)
    res$se     <- msm::deltamethod(resFormula, coef(model), vcov(model))
    res$tvalue <- res$est / res$se
    res$pvalue <- ifelse((type == "t"), ## faire une condition plus intelligente (t and norm)
                         2 * pt(abs(res$tvalue), df.residual(model), lower.tail = FALSE),
                         2 * pnorm(abs(res$tvalue), lower.tail = FALSE))

    # renvoyer les bons noms : par exempel:
    # Estimate Std. Error   z value        Pr(z)
    # ou
    # Estimate Std. Error t value Pr(>|t|)
    res
}


#' @export
extraCoef <- function(model) {
    #myfunc <- function(v1) {
    #    deparse(substitute(v1))
    #}
    #
    #myfunc(foo)
    #[1] "foo"

    # Retrieve result matrice (to which we will add extra estimates)
    coefs   <- coef(stats4::summary(model)) # NON pas normal de devoir mettre bbmle ici, pas besoin de le faire dans estimateCoef
    # not need to do that forcoef(model) ... why?

    # pour betabinom uniquement
    extraCoefs = list(
        p     = quote(x1 / (x1 + x2)),
        theta = quote(1  / (x1 + x2)),
        rho   = quote(1  / (x1 + x2 + 1))
    )

    newCoefs <- do.call(rbind, lapply(extraCoefs, function(i) {
        unlist(estimateCoef(model, i, type = "norm")) # faire plus intelligent ici avec le "norm"
    }))
    rbind(coefs, newCoefs)
}

