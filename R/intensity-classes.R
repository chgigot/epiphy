#------------------------------------------------------------------------------#
#' @include epiphy.R
#' @include utils.R
#------------------------------------------------------------------------------#
NULL

#==============================================================================#
# Definition and basic toolbox for "mapping" objects
#==============================================================================3

#------------------------------------------------------------------------------#
#' Construct data mappings.
#'
#' Data mappings describe how variables in the data are mapped to standard names
#' used throughout \code{epiphy}.
#'
#' Standard names are \code{x}, \code{y} and \code{z} for the three spatial
#' dimensions, and \code{t} for the time. \code{r} corresponds to the records
#' of (disease) intensity, and \code{n}, the number of individuals in a sampling
#' unit (if applicable).
#'
#' \code{mapping()} works with expressions, and \code{mapping_()}, with a vector
#' of characters.
#'
#' @param ... One or more unquoted expressions separated by commas.
#' @param x Vector of one or more character strings.
#' @param data An \code{intensity} object.
#' @param mapping A \code{mapping} object.
#' @param keep_only_std Keep only standard variables.
#'
#' @seealso \code{\link{mapped_var}}
#'
#' @examples
#' mapping(x = col1, y = col2)
#' mapping_(c("x = col1", "y = col2"))
#'
#' @export
#------------------------------------------------------------------------------#
mapping <- function(...) {
    map <- as.list(match.call()[-1])
    map <- lapply(map, function(x) ifelse(is.name(x), x, as.name(x)))
    map <- map[sort(names(map))] # To reorder the names in a standad way.
    structure(map, class = "mapping")
}

#------------------------------------------------------------------------------#
#' @rdname mapping
#' @export
#------------------------------------------------------------------------------#
mapping_ <- function(x) {
    trim       <- function(x) gsub("^[[:space:]]+|[[:space:]]+$", "", x)
    splitted   <- lapply(strsplit(x, "="), trim)
    map        <- lapply(splitted, tail, n = 1L)
    names(map) <- vapply(splitted, head, n = 1L, FUN.VALUE = character(1))
    map        <- lapply(map, as.name)
    map        <- map[sort(names(map))] # To reorder the names in a standad way.
    structure(map, class = "mapping")
}

#------------------------------------------------------------------------------#
#' @rdname mapping
#' @export
#------------------------------------------------------------------------------#
remap <- function(data, mapping, keep_only_std = TRUE) {
    stopifnot(any(class(data) == "intensity"))
    init_intensity(data$data, mapping, keep_only_std, type = is(data))
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.mapping <- function(x, ...) {
    values <- vapply(x, deparse, character(1))
    bullets <- paste0("* ", names(x), " -> ", values, "\n")
    cat(bullets, sep = "")
    invisible(x)
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
str.mapping <- function(object, ...) utils::str(unclass(object), ...)

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
"[.mapping" <- function(x, i, ...) structure(unclass(x)[i], class = "mapping")

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
c.mapping <- function(..., recursive = FALSE) { # To catch and "jail" recursive arg
    structure(c(unlist(lapply(list(...), unclass))), class = "mapping")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
as.character.mapping <- function(x, ...) {
    char <- as.character(unclass(x))
    names(char) <- names(x)
    char
}

#------------------------------------------------------------------------------#
#' @keywords internal
#------------------------------------------------------------------------------#
map_data <- function(object) {
    as.data.frame(lapply(object$mapping, eval, envir = object$data,
                         enclos = NULL), stringsAsFactors = FALSE)
}

#------------------------------------------------------------------------------#
#' @keywords internal
#------------------------------------------------------------------------------#
unmap_data <- function(mapped_data, source_object) {
    flat_mapping <- as.character(source_object$mapping)
    rev_mapping  <- mapping_(paste0(flat_mapping, "=", names(flat_mapping)))
    data <- as.data.frame(lapply(rev_mapping, eval, envir = mapped_data,
                                 enclos = NULL), stringsAsFactors = FALSE)
    structure(list(data    = data,
                   mapping = source_object$mapping,
                   struct  = source_object$struct),
              class = class(source_object))
}

#==============================================================================#
#  Definition of "intensity" objects
# (i.e. "count", "incidence", "severity" objects)
#==============================================================================#

#------------------------------------------------------------------------------#
# Validity of intensity objects
#------------------------------------------------------------------------------#
valid_intensity <- function(object) {

    ## intensity
    if (!any(class(object) == "intensity"))
        stop("Must be an 'intensity' object.")
    if (nrow(object$obs) == 0)
        stop("There must be at least one observation.")
    #if (anyDuplicated(data.frame(object@space, object@time)) != 0)
    #    stop("Data set with duplicated records.")
    # BESOIN DE PRENDRE EN COMPTE LES LABELS / DE MEME POUR LA FUNCTION PRINT, ca doit apparaitre

    ## count
    if (any(class(object) == "count")) {
        if (ncol(object$obs) != 1)
            stop("Each record must have only one observation ('i').")
        if (!all(names(object$obs) %in% "i"))
            stop("Name of observation column must be 'i'.")
        if (!all(is.wholenumber(object$obs)))
            stop("Observation values must be integers.")
        if (!all(object$obs[["i"]] >= 0))
            stop("Observation values must be >= 0")
    }

    ## incidence
    if (any(class(object) == "incidence")) {
        if (ncol(object$obs) != 2)
            stop("Each record must have two observations ('i' and 'n').")
        if (!all(names(object$obs) %in% c("i", "n")))
            stop("Names of observation columns must be 'i' and 'n'.")
        if (!all(is.wholenumber(object$obs)))
            stop("Observation values must be integers.")
    }

    ## severity
    if (any(class(object) == "severity")) {
        if (ncol(object$obs) != 1)
            stop("Each record must have only one observation ('i').")
        if (!all(names(object$obs) %in% "i"))
            stop("Name of observation column must be 'i'.")
        if (!all((object$obs >= 0) & (object$obs <= 1)))
            stop("Observation values must between 0 and 1.")
    }

    TRUE
}

#------------------------------------------------------------------------------#
# Initial checking and building of intensity object
#------------------------------------------------------------------------------#
init_intensity <- function(data, mapping, keep_only_std, type) {

    # TODO: keep_only_std do not work
    # TODO: code to be cleaned a bit.

    std_names <- list(space = c("x", "y", "z"), # up to 3 dim for space
                      time = "t",               # up to 1 dim for time
                      obs = c("i", "n"))        # up to 2 "dim" for observations:
                                                # - r = records
                                                # - n = number of individuals

    #--------------------------------------------------------------------------#
    # Initial checks and attribution of a given code for each class (and struct?)
    # * type
    if (missing(type)) stop("Missing 'type'.")
    # * data
    if (missing(data)) stop("Missing 'data'.")
    if (!(is.data.frame(data))) stop("'data' must be a data frame.")
    # Subsequent verifications should be done in 'validity' function for each
    # object.
    data_header   <- colnames(data) ## Redondant avec plus haut
    # * mapping
    if (missing(mapping)) {
        data_std_names <- data_header[data_header %in% unlist(std_names)]
        mapping <- mapping_(paste0(data_std_names, "=", data_std_names))
        # checkings!!!
    } else {
        if (class(mapping) != "mapping") stop("'mapping' must be a mapping object.")
        names_mapping <- names(mapping) ## Redondant avec plus bas
        if (!all(i_std <- names_mapping %in% unlist(std_names)) && keep_only_std) { # TODO: What? keep_only_std here?
            #warning("Dropping unrelevant names in mapping.") # NOT CLEAR
            mapping <- mapping[i_std]
            # checkings!!!!
        }
        # Then, Check if there are some standard names in colnames that need to be
        # auto-mapped:
        # TODO: faire en sorte que si
        mapped_header <- as.character(mapping)
        unmapped_data_header <- data_header[!(data_header %in% mapped_header)]
        unmapped_data_std_names <- unmapped_data_header[
            (unmapped_data_header %in% unlist(std_names)) & !(unmapped_data_header %in% names_mapping)
            # ^ if it is a standard name ...                ^ and if it is not already mapped
        ]
        if (length(unmapped_data_std_names) > 0) {
            extra_mapping <- mapping_(paste0(unmapped_data_std_names, "=", unmapped_data_std_names))
            mapping <- c(mapping, extra_mapping)
        }
    }

    mapping_names <- names(mapping) ## Redondant avec plus haut
    struct <- lapply(std_names, function(type) {
        mapping_names[mapping_names %in% type]
    })

    # TODO: NEW: to test.
    if (!keep_only_std) {
        unmapped_data_non_std_names <- data_header[!(data_header %in% mapping)]
        extra_mapping <- mapping_(paste0(unmapped_data_non_std_names, "=", unmapped_data_non_std_names))
        mapping <- c(extra_mapping, mapping) ## labels at the beginning (convention)
        struct  <- c(list(label = unmapped_data_non_std_names), struct)  ## labels at the beginning (convention)
    }

    #--------------------------------------------------------------------------#
    # Return an "intensity" object
    structure(list(data    = data,
                   mapping = mapping,
                   struct  = struct),
              class = c(type, "intensity"))
}

#------------------------------------------------------------------------------#
#' Construct count, incidence and severity objects.
#'
#' \code{count()}, \code{incidence()} and \code{severity()} create eponym
#' objects. All of these classes inherit from the base class \code{intensity}.
#' The choice of the class depends on the nature of the data set.
#'
#' \code{incidence} reads disease incidence data from a data frame and return an
#' incidence object. All of these classes inherit from \code{intensity} class.
#' \itemize{
#'     \item count: Each sampling unit contains from 0 to theoreticaly an infinity of data.
#'     Number are positive integers.
#'     \item incidence: Each sampling unit contains an number of diseased plants,
#'     ranging from 0 to \code{n} which is the total amount of plants per sampling
#'     unit.
#'     \item severity: Each sampling unit contain a percentage of disease, a positive
#'     real number ranging from 0.0 to 1.0.
#' }
#'
#' Class intensity and inherited classes
#'
#' All the classes recording disease intensity measurements inherit from this
#' class. The class \code{intensity} is virtual which means that no object of a
#' class \code{intensity} can be constructed. This class only describes common
#' features of all the different disease intensity measurements implemented in
#' this package (\code{\link{count}}, \code{\link{incidence}} and
#' \code{\link{severity}}). You should call one of these inherited classes
#' instead, depending on the nature of your data.
#'
#' By convention, the first columns of the different data frames of each slots
#' have names, but the spatial, temporal or even disease information do not need
#' to fit to these conventions or may be less straightforward and need more
#' columns to record correctly all the information. In such unusual situations,
#' the automatic options of the analysis tools would need to be overridden to be
#' able to work in the desired way.
#'
#' The differences between the different inherited classes regard only the
#' \code{obs} slot. In the case of \code{\link{count}}, the data expected for
#' each record are positive integers (N+). For \code{\link{incidence}}, the data
#' sets are supposed to be two information set per records, the number of
#' diseased unit per sampling unit (r) and the total number of units per
#' sampling unit (n). Note that in its current implementation, n is supposed to
#' be the same for a whole data set. Unequal sampling units are not implemented
#' yet. Finally, for \code{\link{severity}}, r is positive real ranging from 0
#' to 1 and depecting a percentage.
#'
#' space A data frame containing only spatial information. Each row
#'   corresponds to a sampling unit. By convention, the first 3 columns are
#'   names \code{x}, \code{y}, \code{z}.
#'
#' time A data frame containing temporal information. By convention, the
#'   first column is named \code{t}.
#'
#' obs A data frame containing disease observations themselves. The name
#'   of the columns may differ between the sub-class chosed to record the data.
#'
#' Note that it is possible to create a "severity" object but no statistical
#' tools are currently implemented to deal with such an object.
#'
#' An \code{intensity} object contains at very least the "pure" intensity
#' records (column \code{r}) which is a so-called observational variable.
#' Another observational variable, the number of individuals in a sampling unit
#' (\code{n}), is present in the case of a \code{incidence} object. Very often
#' in addition to observational variables, there are spatial (columns \code{x},
#' \code{y} and/or \code{z}) and/or temporal (column \code{t}) variables.
#'
#' Note that the \code{severity} class and the \code{z} variable (the 3rd
#' spatial dimension) are implemented but no statistical methods use them at
#' this point.
#'
#' @param data A data frame. Each line corresponds to a record (or case, or
#'     entry).
#' @param mapping A \code{mapping} object, created with \code{mapping()} or
#'     \code{mapping_()} functions. ... A vector with all the corresponding variables. The different
#'   elements can be named (names of the elements) of the data frame in the
#'   incidence object), or unamed. In the latter case, elements must be
#'   correctly ordered, i.e. x, y, z, t, r and then n. If variables in NULL,
#'   then only the 6 first ... will be take into account in the following (1, 2,
#'   ...), i.e. the id of the value. All the 'parameters' need to be specified.
#' @param keep_only_std Are only standard names kept when proceeding to mapping?
#'   Setting \code{keep_only_std} to TRUE may be useful for subsequent data splitting
#'   using extra labels.
#'
#' @return An \code{incidence} object.
#'
#' When printed, difference information are available:
#'
#' \itemize{
#'     \item The number of sampling units.
#'     \item The time.
#'     \item Is it georeferenced (TRUE/FALSE)
#'     \item Are there any NA data (TRUE/FALSE)
#'     \item Is it a complet array (TRUE/FALSE)? A complete array means that all the recorded values allow to
#'     display an array (even if some data are not available), but this was explicitelly specified. To
#'     complete a dataset, just use \code{complete(data)}. You can also remove NA, which is necessary to use
#'     some analysis technics, using \code{replaceNA(data)} or \code{replace.na(data)}. Note that using both
#'     commands will results in modifying the original data sets which will be specified.
#' }
#'
#'
#' @examples
#' ## Create intensity objects
#' # Implicite call: The variable mapping does not need to be specified if the
#' # column names of the input data frame follow the default names.
#' colnames(tomato_tswv$field_1929) # Returns c("x", "y", "t", "i", "n")
#' my_incidence_1 <- incidence(tomato_tswv$field_1929)
#' my_incidence_1
#' my_incidence_2 <- incidence(tomato_tswv$field_1929,
#'                             mapping(x = x, y = y, t = t, i = i, n = n))
#' identical(my_incidence_1, my_incidence_2)
#'
#' # Explicite call: Otherwise, the variable mapping need to be specified, at
#' # least for column names that do not correspond to default names.
#' colnames(aphids) # Returns c("xm", "ym", "i")
#' my_count_1 <- count(aphids, mapping(x = xm, y = ym, i = i))
#' my_count_1
#' # We can drop the "i = i" in the mapping.
#' my_count_2 <- count(aphids, mapping(x = xm, y = ym))
#' identical(my_count_1, my_count_2)
#'
#' # It is possible to change the variable mapping after the creation of an
#' # intensity object:
#' another_incidence <- incidence(hop_viruses$HpLV)
#' another_incidence
#' remap(another_incidence, mapping(x = xm, y = ym))
#'
#' ## Plotting data
#' plot(my_incidence_1) # Same as: plot(my_incidence_1, type = "spatial")
#' plot(my_incidence_1, type = "temporal")
#'
#' plot(my_count_1, tile = FALSE, size = 5)
#' plot(my_count_1, type = "temporal") # Not possible: there is only 1 date.
#'
#' # Using grayscale:
#' plot(my_count_1, grayscale = TRUE)
#' plot(my_count_1, grayscale = TRUE, tile = FALSE, size = 5)
#'
#' @name intensity
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' @rdname intensity
#' @aliases count_data
#' @export
#------------------------------------------------------------------------------#
count <- function(data, mapping, keep_only_std = TRUE) {
    init_intensity(data, mapping, keep_only_std, type = "count")
}

#------------------------------------------------------------------------------#
#' @rdname intensity
#' @aliases incidence_data
#' @export
#------------------------------------------------------------------------------#
incidence <- function(data, mapping, keep_only_std = TRUE) {
    init_intensity(data, mapping, keep_only_std, type = "incidence")
}

#------------------------------------------------------------------------------#
#' @rdname intensity
#' @aliases severity_data
#' @export
#------------------------------------------------------------------------------#
severity <- function(data, mapping, keep_only_std = TRUE) {
    init_intensity(data, mapping, keep_only_std, type = "severity")
}

# The three following function (*_data) are alternative way to create
# intensity objects. They just are aliases.

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
count_data <- function(data, mapping, keep_only_std = TRUE) {
    count(data, mapping, keep_only_std)
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
incidence_data <- function(data, mapping, keep_only_std = TRUE) {
    incidence(data, mapping, keep_only_std)
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
severity_data <- function(data, mapping, keep_only_std = TRUE) {
    severity(data, mapping, keep_only_std)
}

#==============================================================================#
# Basic manipulations of "intensity" objects
#==============================================================================#

#------------------------------------------------------------------------------#
#' Test if an object is of class \code{intensity} or one of its subclasses.
#'
#' Test if an object is of class \code{intensity} or one of its subclasses
#' (i.e. \code{count}, \code{incidence} or \code{severity}).
#'
#' @param x An object.
#'
#' @export
#------------------------------------------------------------------------------#
is.intensity <- function(x) return(is(x, "intensity"))
#------------------------------------------------------------------------------#
#' @rdname is.intensity
#' @export
#------------------------------------------------------------------------------#
is.count     <- function(x) return(is(x, "count"))
#------------------------------------------------------------------------------#
#' @rdname is.intensity
#' @export
#------------------------------------------------------------------------------#
is.incidence <- function(x) return(is(x, "incidence"))
#------------------------------------------------------------------------------#
#' @rdname is.intensity
#' @export
#------------------------------------------------------------------------------#
is.severity  <- function(x) return(is(x, "severity"))

#------------------------------------------------------------------------------#
# Text outputs for "intensity" objects (print, summary)
#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.intensity <- function(x, ...) {
    # TODO: Clean this function (no variable with name 'test','a', ...).
    cat("# A mapped object: ", class(x)[1], " class\n", sep = "")
    len <- lengths(x$struct)
    cat("# dim: ", len[1], " ", names(x$struct)[1], ", ", len[2], " ", names(x$struct)[2],
        ", ", len[3], " ", names(x$struct)[3], "\n", sep = "")
    heading      <- colnames(x$data)
    mapping_from <- names(x$mapping)
    mapping_to   <- unname(as.character(x$mapping))
    idx_col      <- sapply(heading, function(col) {
        x <- which(mapping_to == col)
        ifelse(length(x) == 1, x, NA)
    })
    test <- as.character(mapping_from[idx_col])
    test[!is.na(test)] <- paste0("[", test[!is.na(test)], "]")
    test[is.na(test)]  <- "."
    test <- setNames(test, heading)
    nrow       <- nrow(x$data)
    nrow_max   <- ifelse(nrow > 6, 6, nrow)
    chr_data   <- x$data[1:nrow_max, ]
    chr_data[] <- lapply(chr_data, as.character) # To avoid problems with factor (labels) column: niveau de facteur incorrect, NAs générés
    a <- data.frame(rbind(names(test), chr_data), row.names = c("", 1:nrow_max))
    colnames(a) <- test
    print(a)
    cat("# ... with ", nrow - nrow_max, " more records (rows)\n", sep = "")
    invisible(x)
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
summary.intensity <- function(object, ...) summary(object$data)

#------------------------------------------------------------------------------#
#' Coerce to a data frame.
#'
#' Functions to coerce an \code{intensity} object to a data frame.
#'
#' @param x An \code{intensity} object.
#' @inheritParams base::as.data.frame
#'
#' @return A data frame.
#'
#' @examples
#' my_data <- incidence(tomato_tswv$field_1929)
#' head(as.data.frame(my_data))
#'
#' @method as.data.frame intensity
#' @export
#------------------------------------------------------------------------------#
as.data.frame.intensity <- function(x, row.names = NULL, optional = FALSE, ...,
                                    stringsAsFactors = default.stringsAsFactors()) {
    # To keep the standard behavior of as.data.frame(), we need to coerce
    # the input data frame to a "simple" list.
    as.data.frame(as.list(x$data), row.names = row.names, optional = optional,
                  ..., stringsAsFactors = stringsAsFactors)
}

#------------------------------------------------------------------------------#
#' Existing variable mappings.
#'
#' Get or set existing variable mappings.
#'
#' @param x An \code{intensity} object.
#' @param value A \code{mapping} object.
#' @param keep Logical. Do we keep any previous mapped variables that are not
#'     redifined in the \code{mapping} object?
#'
#' @seealso \code{\link{mapping}}
#'
#' @examples
#' my_data <- count(aphids)
#' my_data
#... TODO
# x, y, X, Y
#' mapped_var(my_data)
#' mapped_var(my_data) <- mapping(x = X, y = Y)
#' mapped_var(my_data)
#' mapped_var(my_data) <- mapping(x = x, r = r, keep = FALSE)
#' mapped_var(my_data)
#'
#' @export
#------------------------------------------------------------------------------#
mapped_var <- function(x) {
    stopifnot(is.intensity(x))
    x$mapping
}

#------------------------------------------------------------------------------#
#' @rdname mapped_var
#' @export
#------------------------------------------------------------------------------#
"mapped_var<-" <- function(x, keep = TRUE, value) {
    stopifnot(is.intensity(x))
    if (keep) {
        idx_to_keep <- !(names(x$mapping) %in% names(value))
        x$mapping <- structure(c(x$mapping[idx_to_keep], value),
                               class = "mapping")
    } else {
        x$mapping <- value
    }
    x
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
dim.intensity <- function(x) lengths(x$struct)





















#==============================================================================#
# Advanced manipulations of "intensity" objects
#==============================================================================#

#------------------------------------------------------------------------------#
#' Regroup observational data into even clumps of individuals.
#'
#' This function provides a easy way to regroup recorded data into groups of
#' same number of individuals.
#'
#' @param object An \code{intensity} object.
#' @param unit_size Size of a group unit. It must be a named vector, with names
#'     corresponding to non-observational variables (i.e. space and time
#'     variables). If the size of a variable in the data set is not a multiple
#'     of the provided value in \code{unit_size}, some sampling units (the last
#'     ones) will be dropped so that clumps of individuals remain even
#'     throughout the data set.
# TODO: @param group_by Not yet implemented.
#' @param fun Function used to group observational data together.
#' @param inclusive_unspecified Not yet implemented. Do unspecified mapped
#'     variables (different from i and n) need to be included into the bigger
#'     possible sampling unit (TRUE) or splited into as many sampling units as
#'     possible (FALSE, default)?
#' @param ... Additional arguments to be passed to \code{fun}.
#'
#' @examples
#' my_incidence <- incidence(tomato_tswv$field_1929)
#' plot(my_incidence, type = "all")
#'
#' # Different spatial size units:
#' my_incidence_clumped_1 <- clump(my_incidence, unit_size = c(x = 3, y = 3))
#' plot(my_incidence_clumped_1, type = "all")
#'
#' my_incidence_clumped_2 <- clump(my_incidence, unit_size = c(x = 4, y = 5))
#' plot(my_incidence_clumped_2, type = "all")
#'
#' # To get mean disease incidence for each plant over the 3 scoring dates:
#' my_incidence_clumped_3 <- clump(my_incidence, unit_size = c(t = 3), fun = mean)
#' plot(my_incidence_clumped_3)
#'
# Interest of the parameter inclusive_unspecified. TODO: Not yet implemented.
# #my_incidence1 <- clump(my_incidence, unit_size = c(x = 3, y = 3))
# #my_incidence2 <- clump(my_incidence, unit_size = c(x = 3, y = 3, t = 1))
# #identical(my_incidence1, my_incidence2)
#
# #my_incidence3 <- clump(my_incidence, unit_size = c(x = 3, y = 3), inclusive_unspecified = TRUE)
# #my_incidence4 <- clump(my_incidence, unit_size = c(x = 3, y = 3, t = 3))
# #identical(my_incidence3, my_incidence4)
#'
#' @export
#------------------------------------------------------------------------------#
clump <- function(object, ...) UseMethod("clump")

#------------------------------------------------------------------------------#
#' @rdname clump
#' @export
#------------------------------------------------------------------------------#
clump.intensity <- function(object, unit_size, fun = sum,
                            inclusive_unspecified = FALSE, ...) { # TODO: Code inclusive_unspecified

    #--------------------------------------------------------------------------#
    # Initial checks and data preparation
    if (is.null(names(unit_size))) {
        stop("unit_size must be a named vector.")
    }
    mapped_data   <- map_data(object)
    colnames_mapped_data <- colnames(mapped_data)
    obs_names     <- object$struct[["obs"]]
    #non_obs_names <- unname(unlist(object$struct[c("space", "time")]))
    non_obs_names <- colnames_mapped_data[!(colnames_mapped_data %in% obs_names)] # because of if keep_only_std = FALSE
    if (!all(names(unit_size) %in% non_obs_names)) {
        stop(paste0("All unit_size names must exist in mapped data ",
                    "(non-observational types)."))
    }

    #--------------------------------------------------------------------------#
    # Define groups

    # exemple avec incidence(hop_viruses$HpLV, mapping(x=xm,y=ym))
    invisible(lapply(seq_len(length(unit_size)), function(i1) {
        id  <- names(unit_size)[i1]
        val <- mapped_data[[id]]
        val <- (val - min(val)) ## Rescale with 0 as a tmp new base in ordre to Gérer aussi le cas du zér0 !!!
        tmp <- floor(val / unit_size[[id]]) + 1 ## +1 as a new base by convention
        if (length(unique(table(tmp))) > 1) {
            tmp[tmp == max(tmp)] <- NA
        }
        mapped_data[[id]] <<- tmp # mapped_data remains a data frame here.
    }))

    #--------------------------------------------------------------------------#
    # Manage uncomplete non-observational cases
    complete_non_obs_cases <- complete.cases(mapped_data[, non_obs_names])
    if (!all(complete_non_obs_cases)) {
        warning(paste0("To get even clumps of individuals, a total of ",
                       sum(!complete_non_obs_cases),
                       " source sampling units were dropped."))
    }
    mapped_data <- mapped_data[complete_non_obs_cases, ]

    #--------------------------------------------------------------------------#
    # Split the data frame
    f <- mapped_data[, non_obs_names]
    # There is no more NA in f, but they may still remain some NA in
    # observational data (i.e. in r and n? columns)
    split_data <- base::split(mapped_data, f = f, drop = TRUE) #, lex.order = TRUE)... only in most recent R version # Only cosmetic (to get a nice order)

    # Make the real calculation
    clumped_data <- do.call(rbind, lapply(split_data, function(sub_data) {
        lapply(colnames(sub_data), function(colname) {
            if (colname %in% non_obs_names) {
                sub_data[[colname]][1]
            } else {
                # ie. if (colname %in% obs_names)
                fun(sub_data[[colname]], ...)
            }
        })
    }))
    # Rearrange a bit the data
    # Note: We need to do all of that to avoid any conversion to a matrix
    # just in case we have different types of data (e.g. numeric and character).
    clumped_data <- setNames(lapply(seq_len(ncol(clumped_data)),
                                 function(i1) unname(unlist(clumped_data[, i1]))),
                             colnames_mapped_data)
    #--------------------------------------------------------------------------#
    # Return an "intensity" object
    unmap_data(droplevels(clumped_data), source_object = object)
    ## above droplevels just in case...
}

#------------------------------------------------------------------------------#
#' Divide into groups and reassemble.
#'
#' Divide into groups and reassemble.
#'
#' @inheritParams base::split
#' @param by The name(s) of the variable(s) which define(s) the grouping.
#' @param unit_size Size of a group unit. It must be a named vector, with names
#'     corresponding to non-observational variables (i.e. space and time
#'     variables). If the size of a variable in the data set is not a multiple
#'     of the provided value in \code{unit_size}, some sampling units (the last
#'     ones) will be dropped so that clumps of individuals remain even
#'     throughout the data set.
#'
# @examples
# #inc_spl_t <- split(inc_clu, by = "t")
# #inc_spl_tbis <- split(inc_clu, unit_size = c(x = 8, y = 20, t = 1))
# #identical(unname(inc_spl_t), unname(inc_spl_tbis))
#'
#' @export
#------------------------------------------------------------------------------#
split.intensity <- function(x, f, drop = FALSE, ..., by, unit_size) {
    mapped_data <- map_data(x)
    if (!missing(by) && !missing(unit_size)) {
        stop("'by' and 'unit_size' cannot be provided at the same time.")
    }
    if (!missing(by)) {
        if (!missing(f)) stop("f and by cannot be given at the same time.")
        stopifnot(all(by %in% names(x$mapping)))
        f <- lapply(by, function(var) getElement(mapped_data, var))
    }
    if (!missing(unit_size)) {
#==============================================================================#
        ## TODO: BEG tmp
        object <- x
        ## TODO: END tmp
        #--------------------------------------------------------------------------#
        # Initial checks and data preparation
        if (is.null(names(unit_size))) {
            stop("unit_size must be a named vector.")
        }
        mapped_data   <- map_data(object)
        colnames_mapped_data <- colnames(mapped_data)
        obs_names     <- object$struct[["obs"]]
        #non_obs_names <- unname(unlist(object$struct[c("space", "time")]))
        non_obs_names <- colnames_mapped_data[!(colnames_mapped_data %in% obs_names)] # because of if keep_only_std = FALSE
        if (!all(names(unit_size) %in% non_obs_names)) {
            stop(paste0("All unit_size names must exist in mapped data ",
                        "(non-observational types)."))
        }

        #--------------------------------------------------------------------------#
        # Define groups

        # TODO : Gérer aussi le cas du zér0 !!!
        # exemple avec incidence(hop_viruses$HpLV, mapping(x=xm,y=ym))
        invisible(lapply(seq_len(length(unit_size)), function(i1) {
            id  <- names(unit_size)[i1]
            tmp <- ceiling(mapped_data[[id]] / unit_size[[id]])
            if (length(unique(table(tmp))) > 1) {
                tmp[tmp == max(tmp)] <- NA
            }
            ## TODO: BEG NEW
            mapped_data[[paste0("tmp_", id)]] <<- mapped_data[[id]]
            ## TODO: END NEW
            mapped_data[[id]] <<- tmp # mapped_data remains a data frame here.
        }))

        #--------------------------------------------------------------------------#
        # Manage uncomplete non-observational cases
        complete_non_obs_cases <- complete.cases(mapped_data[, non_obs_names])
        if (!all(complete_non_obs_cases)) {
            warning(paste0("To get even clumps of individuals, a total of ",
                           sum(!complete_non_obs_cases),
                           " source sampling units were dropped."))
        }
        mapped_data <- mapped_data[complete_non_obs_cases, ]

        #--------------------------------------------------------------------------#
        # Split the data frame
        f <- mapped_data[, non_obs_names]
        # There is no more NA in f, but they may still remain some NA in
        # observational data (i.e. in r and n? columns)
        #split_data <- base::split(mapped_data, f = f) #, lex.order = TRUE)... only in most recent R version # Only cosmetic (to get a nice order)
#==============================================================================#
    }
    split_data <- base::split(mapped_data, f, drop, ...)
#==============================================================================#
    if (!missing(unit_size)) {
        #--------------------------------------------------------------------------#
        # Rename the variables ## TODO: NEW
        split_data <- lapply(split_data, function(sub_data) {
            colnames_sub_data <- colnames(sub_data)
            tmp_vars <- colnames_sub_data[grep("tmp_", colnames_sub_data)]
            lapply(tmp_vars, function(tmp_var) {
                new_var <- sub("tmp_", "", tmp_var)
                sub_data[[new_var]] <<- sub_data[[tmp_var]]
                sub_data[[tmp_var]]             <<- NULL
            })
            sub_data
        })
    }
#==============================================================================#

    ## below droplevels just in case there were factors that need to be "cleaned" aftre the spliting
    lapply(split_data, function(subx) unmap_data(droplevels(subx), source_object = x))
}

# TODO: unsplit

# TODO: Doc and clean below

#------------------------------------------------------------------------------#
#' To go to higher level in the hierarchy.
#'
#' This function transforms the current numeric vector or \code{intensity} data
#' set into a "simplified black and white image" of this same data set: every
#' value of disease intensity below and above a given threshold is given the
#' value 0 and 1, respectively.
#'
#' By default, everything above 0 is given 1, and 0 stays at 0. \code{threshold}
#' is thus useful to report a whole sampling unit as "healthy" (0), if no
#' diseased individual at all was found within the sampling unit, or "diseased"
#' (1) if at least one diseased individual was found.
#'
#' @param data A numeric vector or an \code{intensity} object.
#' @param value All the intensity values lower or equal to this value  are set
#'     to 0. The other values are set to 1.
#' @param ... Additional arguments to be passed to other methods.
#'
#' @export
#------------------------------------------------------------------------------#
threshold <- function(data, value, ...) UseMethod("threshold")

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
threshold.numeric <- function(data, value = 0, ...) {
    ifelse(data <= value, 0, 1)
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
threshold.intensity <- function(data, value = 0, ...) {
    mapped_data <- map_data(data)
    mapped_data[["i"]] <- ifelse((mapped_data[["i"]] <= value), 0, 1)
    if ("n" %in% colnames(mapped_data)) {
        mapped_data[["n"]] <- 1
    }
    unmap_data(mapped_data, data)
}

#------------------------------------------------------------------------------#
# TODO: to keep below?
is.completeArray <- function(object) {
    dt <- as.data.frame(object, fields = c("space", "time"))
    ref <- expand.grid(lapply(dt, unique))
    #names(dt) <- names(ref) <- make.names(rep("X", ncol(dt)), unique = TRUE)
    res <- merge(dt, ref, all = TRUE)
    if (nrow(res) == nrow(dt)) return(TRUE)
    else return(FALSE)
}

#------------------------------------------------------------------------------#
# Miscellaneous
#------------------------------------------------------------------------------#
# Useful to display in a data.frame de type list-column / column-list
#------------------------------------------------------------------------------#

# Apparently not needed to create setGeneric("as.character")
# Seems to be automatically done when creating setMethdod(...)
# To double-check
# @export
# setGeneric("as.character")
# Actually, it seems to generate a warning:
# Warning message:
# In setup_ns_exports(pkg, export_all) :
#    Objects listed as exports, but not present in namespace: as.character

#------------------------------------------------------------------------------#
#' @method as.character count
#' @export
#------------------------------------------------------------------------------#
as.character.count <- function(x, ...) "<count object>"

#------------------------------------------------------------------------------#
#' @method as.character incidence
#' @export
#------------------------------------------------------------------------------#
as.character.incidence <- function(x, ...) "<incidence object>"

#------------------------------------------------------------------------------#
#' @method as.character severity
#' @export
#------------------------------------------------------------------------------#
as.character.severity <- function(x, ...) "<severity object>"


#==============================================================================#
# Graphical outputs for "intensity" objects
#==============================================================================#

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
plot.intensity <- function(x, y, ..., type = c("spatial", "temporal", "all"),
                           tile = TRUE, pch = 22, legend.position = "right",
                           grayscale = FALSE) {

    type <- match.arg(type)

    #TODO: if (!missing(y)) ...
    #mapping <- attr(x, "mapping")
    mapped_data <- map_data(x)
    label       <- lapply(x$mapping, deparse) # TODO: useful?

    possible_temporal  <- ("t" %in% colnames(mapped_data)) # TODO: More security here
    possible_x_spatial <- ("x" %in% colnames(mapped_data)) # TODO: More security here
    possible_y_spatial <- ("y" %in% colnames(mapped_data)) # TODO: More security here
    possible_spatial   <- possible_x_spatial || possible_y_spatial # TODO: More security here

    max_scale <- ifelse(class(x)[1] == "incidence",
                        max(mapped_data[["n"]]),
                        max(mapped_data[["i"]]))
    fill_breaks <- if (max_scale <= 10) seq(0, max_scale)  else waiver()
    fill_limits <- if (max_scale <= 10) range(fill_breaks) else c(0, max_scale)

    color_high <- ifelse(grayscale, "black", "red")

    # Temporal figure
    if (type %in% c("temporal", "all") && possible_temporal) {
        nt <- length(unique(mapped_data$t))
        ni <- length(unique(mapped_data$i))
        t_breaks <- if(nt <= 10) unique(mapped_data$t) else waiver()
        i_breaks <- if(ni <= 10) unique(mapped_data$i) else waiver()
        gg <- ggplot()
        gg <- gg + geom_jitter(data = mapped_data, mapping = aes(t, i),
                               alpha = 0.2, width = 0.2, height = 0)
        gg <- gg + stat_summary(data = mapped_data, mapping = aes(t, i),
                                fun.y = "mean", geom = "line", color = color_high,
                                linetype = "dashed")
        gg <- gg + stat_summary(data = mapped_data,
                                mapping = aes(t, i, group = t),
                                fun.data = "mean_sdl", fun.args = list(mult = 1),
                                geom = "pointrange", # default
                                color = color_high)
        gg <- gg + scale_x_continuous("Time (t)", breaks = t_breaks)
        gg <- gg + scale_y_continuous(paste0(tocamel(class(x)[1]), " (i)"),
                                      breaks = i_breaks)
        gg <- gg + theme_bw()
        print(gg)
    }

    # Spatial figure
    if (type %in% c("spatial", "all") && possible_spatial) {
        nx <- length(unique(mapped_data$x))
        ny <- length(unique(mapped_data$y))
        x_breaks <- if(nx <= 10) unique(mapped_data$x) else waiver()
        y_breaks <- if(ny <= 10) unique(mapped_data$y) else waiver()
        gg <- ggplot()
        if (tile) {
            gg <- gg + geom_raster(data = mapped_data,
                                   mapping = aes(x, y, fill = i), ...)
            gg <- gg + scale_x_continuous("Spatial dim 1 (x)",
                                          breaks = x_breaks, expand = c(0, 0))
            gg <- gg + scale_y_continuous("Spatial dim 2 (y)",
                                          breaks = y_breaks, expand = c(0, 0))
        } else {
            gg <- gg + geom_point(data = mapped_data,
                                  mapping = aes(x, y, fill = i), pch = pch, ...)
            gg <- gg + scale_x_continuous("Spatial dim 1 (x)",
                                          breaks = x_breaks)
            gg <- gg + scale_y_continuous("Spatial dim 2 (y)",
                                          breaks = y_breaks)
        }
        gg <- gg + scale_fill_gradient(paste0(tocamel(class(x)[1]), " (i)"),
                                       low = "white", high = color_high,
                                       guide = guide_legend(reverse = TRUE),
                                       breaks = fill_breaks,
                                       limits = fill_limits)
        gg <- gg + theme_bw() +
            theme(panel.grid = element_blank()) +
            theme(legend.position = legend.position)
        gg <- gg + coord_fixed() # To get squares in any case
        if (possible_temporal) {
            nt <- length(unique(mapped_data$t))
            gg <- gg + facet_wrap(~ t, labeller = label_both, # TODO: Improve labeller to add "Time" everywhere
                                  # To get a pretty output:
                                  nrow = ifelse(nt <= 3, 1, floor(sqrt(nt))))
        }
        print(gg)
    }

    invisible(NULL)

}




