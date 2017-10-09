#------------------------------------------------------------------------------#
#' @include epiphy.R
#' @include utils.R
#------------------------------------------------------------------------------#
NULL

#==============================================================================#
# Definition and basic toolbox for "mapping" objects
#==============================================================================3

#------------------------------------------------------------------------------#
#' Construct data mappings
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
#'
#' @seealso \code{\link{mapv}}
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
    structure(map, class = "mapping")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.mapping <- function(x, ...) {
    values <- vapply(x, deparse, character(1))
    bullets <- paste0("* ", names(x), " -> ", values, "\n")
    cat(bullets, sep = "")
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
            stop("Each record must have only one observation ('r').")
        if (!all(names(object$obs) %in% "r"))
            stop("Name of observation column must be 'r'.")
        if (!all(is.wholenumber(object$obs)))
            stop("Observation values must be integers.")
        if (!all(object$obs$r >= 0))
            stop("Observation values must be >= 0")
    }

    ## incidence
    if (any(class(object) == "incidence")) {
        if (ncol(object$obs) != 2)
            stop("Each record must have two observations ('r' and 'n').")
        if (!all(names(object$obs) %in% c("r", "n")))
            stop("Names of observation columns must be 'r' and 'n'.")
        if (!all(is.wholenumber(object$obs)))
            stop("Observation values must be integers.")
    }

    ## severity
    if (any(class(object) == "severity")) {
        if (ncol(object$obs) != 1)
            stop("Each record must have only one observation ('r').")
        if (!all(names(object$obs) %in% "r"))
            stop("Name of observation column must be 'r'.")
        if (!all((object$obs >= 0) & (object$obs <= 1)))
            stop("Observation values must between 0 and 1.")
    }

    TRUE
}

#------------------------------------------------------------------------------#
# Initial checking and building of intensity object
#------------------------------------------------------------------------------#
init_intensity <- function(data, mapping, type) {

    std_names <- list(space = c("x", "y", "z"), # up to 3 dim for space
                      time = "t",               # up to 1 dim for time
                      obs = c("r", "n"))        # up to 2 "dim" for observations:
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
    # * mapping
    if (missing(mapping)) {
        data_header     <- colnames(data)
        data_std_names <- data_header[data_header %in% unlist(std_names)]
        mapping <- mapping_(paste0(data_std_names, "=", data_std_names))
        # checkings!!!
    } else {
        if (class(mapping) != "mapping") stop("'mapping' must be a mapping object.")
        names_mapping <- names(mapping)
        if (!all(i_std <- names_mapping %in% unlist(std_names))) {
            #warning("Dropping unrelevant names in mapping.") # NOT CLEAR
            mapping <- mapping[i_std]
            # checkings!!!!
        }
    }

    mapping_names <- names(mapping)
    struct <- lapply(std_names, function(type) {
        mapping_names[mapping_names %in% type]
    })

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
#' incidence object. All of these classes inherit from \code{\link{intensity-class}}.
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
#' @seealso \code{\link{intensity-class}}
#' @examples
#' # Implicite call: The parameter struct does not need to be specified if
#' # the column names of the input data frame respect the convention.
#' colnames(Cochran1936) # Returns c("x", "y", "t", "r", "n")
#' incidence(Cochran1936)
#'
#' # Explicit call: Otherwise, struct must be present:
#' incidence(Cochran1936, c(r = "r", n = "n", t = "t", x = "x", y = "y"))
#' incidence(Cochran1936, c(r = 4, n = 5, t = 3, x = 1, y = 2))
#' incidence(Cochran1936, list(space = 1:2, time = 3, obs = 4:5))
#'
#' ## If a variable is not specified, this means it does not exist in the
#' ## input data frame.
#' subData <- subset(Cochran1936, t == 1,
#'                   select = c("x", "y", "r", "n"))
#' # The two following instructions work:
#' incidence(subData)
#' incidence(subData, c(x = 1, y = 2, r = 4, n = 5))
#'
#' records <- incidence(tomato_tswv$field2)
#' records
#' summary(records)
#' plot(records, type = "map_2d", t = 1)
#'
#' ##--
#' my_raw_data <- tomato_tswv$field_1929
#' names(my_raw_data)
#' # Automatic mapping because the name of ...
#' my_data <- incidence(my_raw_data)
#' # To make it less trivial, let's change the column names:
#' colnames(my_raw_data) <- c("coord1", "coord2", "time", "scoring", "tot_plants")
#' my_mapping <- mapping(x = coord1, y = coord2, t = time, r = scoring, n = tot_plants)
#' my_mapping
#' my_data <- incidence(my_raw_data, my_mapping)
#' my_data
#' summary(my_data)
#' plot(my_data)
#'
#' @name intensity
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' @rdname intensity
#' @export
#------------------------------------------------------------------------------#
count <- function(data, mapping) {
    init_intensity(data, mapping, type = "count")
}

#------------------------------------------------------------------------------#
#' @rdname intensity
#' @export
#------------------------------------------------------------------------------#
incidence <- function(data, mapping) {
    init_intensity(data, mapping, type = "incidence")
}

#------------------------------------------------------------------------------#
#' @rdname intensity
#' @export
#------------------------------------------------------------------------------#
severity <- function(data, mapping) {
    init_intensity(data, mapping, type = "severity")
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
    a <- data.frame(rbind(names(test), x$data[1:6, ]), row.names = c("", 1:6))
    colnames(a) <- test
    print(a)
    cat("# ... with ", nrow(x$data) - 6, " more records (rows)\n")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
summary.intensity <- function(object, ...) summary(object$data)

#------------------------------------------------------------------------------#
#' Coerce to a data frame
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
#' Existing variable mappings
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
#' #... TODO
#' # x, y, X, Y
#' mapv(my_data)
#' mapv(my_data) <- mapping(x = X, y = Y)
#' mapv(my_data)
#' mapv(my_data) <- mapping(x = x, r = r, keep = FALSE)
#' mapv(my_data)
#'
#' @export
#------------------------------------------------------------------------------#
mapv <- function(x) {
    stopifnot(is.intensity(x))
    x$mapping
}

#------------------------------------------------------------------------------#
#' @rdname mapv
#' @export
#------------------------------------------------------------------------------#
"mapv<-" <- function(x, value, keep = TRUE) {
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
#' Regroup observational data into even clumps of individuals
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
#' @param group_by Not yet implemented.
#' @param fun Function used to group observational data together.
#' @param ... Optional arguments to \code{fun}.
#'
#' @examples
#' my_data1 <- incidence(tomato_tswv$field_1929)
#' summary(my_data1)
#' plot(my_data1)
#'
#' my_data2 <- clump(my_data1, unit_size = c(x = 3, y = 3))
#' summary(my_data2)
#' plot(my_data2)
#'
#' my_data3 <- clump(my_data1, unit_size = c(t = 3), fun = mean)
#' summary(my_data3)
#' plot(my_data3)
#'
#' @export
#------------------------------------------------------------------------------#
clump <- function(object, ...) UseMethod("clump")

#------------------------------------------------------------------------------#
#' @rdname clump
#' @export
#------------------------------------------------------------------------------#
clump.intensity <- function(object, unit_size, fun = sum, ...) {

    #--------------------------------------------------------------------------#
    # Initial checks and data preparation
    if (is.null(names(unit_size))) {
        stop("unit_size must be a named vector.")
    }
    non_obs_names <- unname(unlist(object$struct[c("space", "time")]))
    obs_names     <- object$struct[["obs"]]
    mapped_data   <- map_data(object)
    if (!all(names(unit_size) %in% non_obs_names)) {
        stop(paste0("All unit_size names must exist in mapped data ",
                    "(non-observational types)."))
    }

    #--------------------------------------------------------------------------#
    # Define groups
    invisible(lapply(seq_len(length(unit_size)), function(i) {
        id  <- names(unit_size)[i]
        tmp <- ceiling(mapped_data[[id]] / unit_size[[id]])
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
    split_data <- base::split(mapped_data, f = f, lex.order = TRUE)

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
                                 function(i) unname(unlist(clumped_data[, i]))),
                             colnames(mapped_data))
    #--------------------------------------------------------------------------#
    # Return an "intensity" object
    unmap_data(clumped_data, source_object = object)
}

#------------------------------------------------------------------------------#
#' Divide into groups and reassemble
#'
#' TODO
#'
#' @inheritParams base::split
#'
#' @export
#------------------------------------------------------------------------------#
split.intensity <- function(x, f, drop = FALSE, ..., by) {
    if (!missing(by)) {
        if (!missing(f)) stop("'f' and 'by' cannot be given at the same time.")
        stopifnot(all(by %in% names(x$mapping)))
        f <- lapply(by, function(var) getElement(x$data, var))
    }
    res <- split(x$data, f, drop, ...)
    lapply(res, function(subx) {
        init_intensity(subx, mapping = x$mapping, type = class(x)[1L])
    })
}

# TODO: unsplit

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
plot.intensity <- function(x, y, ..., type = c("spatial", "temporal", "both"),
                           tile = TRUE, pch = 22) {

    type <- match.arg(type)

    #if (!missing(y)) ...
    #mapping <- attr(x, "mapping")
    mapped_data <- map_data(x)
    label       <- lapply(x$mapping, deparse) # useful?

    max_fill_scale <- ifelse(class(x)[1] == "incidence",
                             max(mapped_data$n),
                             max(mapped_data$r))

    gg <- list()

    #label_full <- ... with time, space, obs...

    # For incidence:
    #seq_gdt <- scales::seq_gradient_pal(low = "white", high = "red") # seq_gdt is a function here
    #seq_gdt <- seq_gdt(seq(0, 1, length = max_fill_scale + 1))

    # Spatial figure
    if (type %in% c("spatial", "both")) {

        # List of layers
        gg_sub <- list(
            geom_raster(data = mapped_data,
                        mapping = aes(x, y, fill = r), ...),
            geom_point(data = mapped_data,
                       mapping = aes(x, y, fill = r), pch = pch, ...),
            #scale_fill_manual(name = paste0(class(x)[[1]], " (r)"),
            #                  values = seq_gdt,
            #                  guide = guide_legend(reverse = TRUE)),
            scale_fill_gradient(name = paste0(class(x)[1], " (r)"),
                                low = "white", high = "red",
                                guide = guide_legend(reverse = TRUE),
                                #breaks = seq(0, max_fill_scale), # To improve later
                                limits = c(0, max_fill_scale)),
            scale_x_continuous(breaks = 1:max(mapped_data$x), expand = c(0, 0)),
            scale_y_continuous(breaks = 1:max(mapped_data$y), expand = c(0, 0)),
            #theme_bw(),
            theme(panel.grid = element_blank()),
            coord_fixed(), # Pour avoir des carr'es pour sûr
            facet_wrap(~ t, labeller = label_both)#, # to wrap a 1d ribbon of panels into 2d.

            ## Valeurs fixées à généraliser ci-après !!!!!!!!
            #scale_colour_gradient(name = "Disease\nintensity",
            #                      low = "white", high = "red", breaks = seq(0, maxScale),
            #                      limits = c(0, maxScale), guide = "legend"),
            # Faire un conditionel ici, de la forme:
            # si pas "reverse x", alors:

            # sinon
            #scale_x_reverse(breaks = xScale, expand = c(0,0)),
            #

            #scale_y_discrete(breaks = yScale, expand = c(0,0)),

            #                         expand = c(0,0)),
            ##
            ## expand = c(0,0) : - OK avec geom_tile
            ##                   - NON avec geom_point
            ##

            # theme(panel.grid = element_blank()),

            #coord_flip(), # Le faire en conditionnel, avec formula : x ~ y (default), ou y ~ x
            #--coord_fixed(), # Pour avoir des carr'es pour sûr
            #--facet_wrap(~ t, nrow = nRibbon)#, # to wrap a 1d ribbon of panels into 2d.
            ### to generaliser pour nrow et ncol
        )

        # Option management
        to_select <- !vector("logical", length = length(gg_sub)) # To create a vector of true
        if (tile) to_select[2] <- FALSE
        else      to_select[1] <- FALSE

        # Store figure information
        gg[[length(gg) + 1]] <- ggplot() + gg_sub[to_select]

    }

    # Temporal figure
    if (type %in% c("temporal", "both")) {

        # List of layers
        gg_sub <- list(
            geom_jitter(data = mapped_data,
                        mapping = aes(t, r), alpha = 0.2,
                        width = 0.2, height = 0),
            stat_summary(data = mapped_data,
                         mapping = aes(t, r),
                         fun.y = "mean", geom = "line", color = "red",
                         linetype = "dashed"),
            stat_summary(data = mapped_data,
                         mapping = aes(t, r, group = t),
                         fun.data = "mean_sdl", fun.args = list(mult = 1),
                         geom = "pointrange", # default
                         color = "red"),
            scale_x_continuous(breaks = seq(0: max(mapped_data$t))),
            scale_y_continuous(breaks = seq(0, max(mapped_data$r))),
            # scale_y_continuous("score (1-9 scale)", breaks=0:9), #### à Généralsier !!!!7        ######### De plus attention, warning quand on combine les graphs
            #### Scale for 'y' is already present. Adding another scale for 'y', which will replace the existing scale.
            #g <- g + scale_x_continuous("date", breaks=seq(1,i))
            expand_limits(y = range(mapped_data$r)),
            theme_bw()
        )

        #scale_y_continuous("score (1-9 scale)", breaks=0:9), #### à Généralsier !!!!
        ######### De plus attention, warning quand on combine les graphs
        #### Scale for 'y' is already present. Adding another scale for 'y', which will replace the existing scale.
        #### g <- g + scale_x_continuous("date", breaks=seq(1,i))
        #geom_point(data = data, aes(x = t, y = r),
        #           alpha = 0.2, size = 3, colour = col,
        #           position = position_jitter(w = 0.2, h = 0.2)),
        #geom_line(data = recap, aes(x = t, y = meanD),
        #          size = 1.5, linetype = 2, colour = col),
        #geom_point(data = recap, aes(x = t, y = meanD),
        #           size = 6, shape = 15, colour = col)

        # Option management
        to_select <- !vector("logical", length = length(gg_sub)) # To create a vector of true
        # Nothing at this point.

        # Store figure information
        gg[[length(gg) + 1]] <- ggplot() + gg_sub[to_select]
    }

    # Dispaly figures (nice to do that at the end to potentialy expand possibilities with + theme... e.g.)
    if (length(gg) == 1) gg[[1]] # Return a gg object
    else                 rev(gg)      # Return a list of gg objects with space fig at the end

}




