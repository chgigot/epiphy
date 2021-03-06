% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intensity-classes.R
\name{split.intensity}
\alias{split.intensity}
\title{Divide into groups and reassemble.}
\usage{
\method{split}{intensity}(x, f, drop = FALSE, ..., by, unit_size)
}
\arguments{
\item{x}{vector or data frame containing values to be divided into groups.}

\item{f}{a \sQuote{factor} in the sense that \code{\link{as.factor}(f)}
    defines the grouping, or a list of such factors in which case their
    interaction is used for the grouping.}

\item{drop}{logical indicating if levels that do not occur should be dropped
    (if \code{f} is a \code{factor} or a list).}

\item{...}{further potential arguments passed to methods.}

\item{by}{The name(s) of the variable(s) which define(s) the grouping.}

\item{unit_size}{Size of a group unit. It must be a named vector, with names
corresponding to non-observational variables (i.e. space and time
variables). If the size of a variable in the data set is not a multiple
of the provided value in \code{unit_size}, some sampling units (the last
ones) will be dropped so that clumps of individuals remain even
throughout the data set.}
}
\description{
Divide into groups and reassemble.
}
