#------------------------------------------------------------------------------#
#' @include epiphy.R
#' @include utils.R
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
mapping <- function(...) {
    map <- as.list(match.call()[-1])
    map <- lapply(map, function(x) ifelse(is.name(x), x, as.name(x)))
    structure(map, class = "mapping")
}

#------------------------------------------------------------------------------#
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

map_data <- function(x) {
    stopifnot("intensity" %in% class(x))
    lapply(x$mapping, eval, envir = x$data, enclos = NULL)
}


#==============================================================================#
# Validity of intensity objects
#==============================================================================#
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

#==============================================================================#
# Initial checking and building of intensity object
#==============================================================================#
init_intensity <- function(data, mapping, type) {

    std_header <- c("x", "y", "z", # up to 3 dim for space
                    "t",           # up to 1 dim for time
                    "r", "n")      # up to 2 "dim" for observations:
                                   # - r = record
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
        data_std_header <- data_header[data_header %in% std_header]
        mapping <- mapping_(paste0(data_std_header, "=", data_std_header))
        # checkings!!!
    } else {
        if (class(mapping) != "mapping") stop("'mapping' must be a mapping object.")
        names_mapping <- names(mapping)
        if (!all(i_std <- names_mapping %in% std_header)) {
            #warning("Dropping unrelevant names in mapping.") # NOT CLEAR
            mapping <- mapping[i_std]
            # checkings!!!!
        }
    }
    # build obj
    #std_data <- std[std$header %in% names(mapping), ]

    # slice <- function(data, std_data, subdf) {
    #     subset(data, select = std_data[std_data$slot == subdf, "header"])
    # }
    #
    # res <- lapply(c("label", "space", "time", "obs"), function(subdf) {
    #     if (subdf == "label") {
    #         `%out%` <- function(x, table) match(x, table, nomatch = 0L) == 0L
    #         res <- subset(data, select = colnames(data)[colnames(data) %out% std_data$header])
    #     } else {
    #         res <- slice(data, std_data, subdf)
    #     }
    #     # In case of empty data frame
    #     if (any(dim(res) == 0)) res <- NULL
    #     res
    # })

    object <- structure(list(data    = data,
                             mapping = mapping),
                        class = c(type, "intensity"))

    #####valid_intensity(object)

    object
    # my_data <- intensity(mtcars, mapping(x = mpg, y = disp))


}

plot.intensity <- function(x, y, ...) {
    #if (!missing(y)) ...
    #mapping <- attr(x, "mapping")
    mapped_data   <- lapply(x$mapping, eval, envir = x$data, enclos = NULL)
    label <- lapply(x$mapping, deparse)
    plot(x    = mapped_data$x,
         y    = mapped_data$y,
         xlab = label$x,
         ylab = label$y)
}

# plot(my_data)


initIntensity <- function(data, struct = NULL, type) { # Find something else that data???????

    ### Split and reshape input data to fit classes input form
    convention <- data.frame(slots   = c(rep("space", 3), "time", rep("obs", 2)), # Find something else that slot???
                             headers = c("x", "y", "z", "t", "r", "n"),
                             stringsAsFactors = FALSE)

    #--------------------------------------------------------------------------#
    # Initial checks and attribution a a given code for each class and struct

    # * type
    if (missing(type)) stop("type must be specified")

    # * data
    if (missing(data)) stop("No data to work with.")
    if (!(is.data.frame(data))) stop("data must be a data frame.")
    #### Subsequent verification will be done in validity for each obj

    # * struct
    if (is.null(struct)) {
        # * struct attribute in data, useful for transformation and retro-transforamtion between data frame and incidence
        if (!is.null(attr(data, "struct"))) {
            stop("En cours...")
        } else {
            dataNamesNoNa <- names(data)[!is.na(names(data))]
            #if (!all(dataNamesNoNa %in% convention$headers))
            #    stop(paste0("If struct is not specified, data headers must ",
            #                "respect the convention."))
            if (length(dataNamesNoNa) != length(unique(dataNamesNoNa)))
                stop(paste0("Data headers must be",
                            "different from each others (except with NA).")) # Faire en sorte que si NA, nom retiré... ou alors retirer cette option

            struct <- names(data)
            names(struct) <- struct
            structType <- "structCharChar"
        }

    } else { # If struct != NULL
        if(!is.vector(struct)) # Note: lists are also considered as vectors in R.
            stop("When specified, struct must be a regular vector or a list.")

        if (is.list(struct)) { # If struct is a list
            flattenstruct <- unlist(struct)
            if (!all(is.wholenumber(flattenstruct)))
                stop(paste0("If struct is a list, all elements must contains",
                            "integer corresponding to data columns indices."))
            if (any(flattenstruct < 0) || any(flattenstruct > ncol(data)))
                stop("BAD indices")
            if (any(is.na(flattenstruct)))
                stop("NA are not allowded in struct, when struct is a list.")
            if (anyDuplicated(flattenstruct))
                stop("At least one index was present twice in struct.")

            structNames <- names(struct)
            if (!all(structNames %in% c(unique(convention$slots), "labels"))) # PAS tres propre du tout ici
                stop(paste0("struct names must be fit with the slot convention:",
                            convention$slots, collapse = ", ")) ## PAS PRPORE ICI, à vérifier
            ##### Ajouter 0 comme indince si d2 par exemple n'est pas pr'esent

            structType <- "structListInt"

        } else if (is.character(struct)) {# struct as a vector of caracters
            # ici, les NA sont autoriser (seul cas !!!!)
            ####if (length(struct) != ncol(data)) stop("BAD1") ... bof bof, car on peut mettre
            # une plus grosse data frame et ne selectionner que qq colonnes
            if (!all(names(struct) %in% convention$headers)) stop("BAD2")
            structNoNa <- struct[!is.na(struct)]
            if (length(structNoNa) != length(unique(structNoNa))) stop("BAD3")
            structType <- "structCharChar"

        } else if (all(is.wholenumber(struct))) { # struct as a vector of integers
            if (any(struct < 0) || any(struct > ncol(data))) stop("BAD indices")
            structNames <- names(struct)
            if (is.null(structNames)) stop("BAD5")
            if (any(is.na(structNames))) stop("BAD5")
            if (!all(structNames %in% convention$headers)) stop("BAD6")
            if (length(structNames) != length(unique(structNames))) stop("BAD7")
            ##### Ajouter 0 comme indince si d2 par exemple n'est pas pr'esent
            structType <- "structCharInt"

        } else {
            stop("struct does not respect conventions.")
        }
    }

    #--------------------------------------------------------------------------#
    # Decoupe proprement !!! en fonction du code (et de la class) ?
    as.indices <- function(slotName, dataHeaders, struct) {
        # Function to extract column indices in the proper order for only struct
        # types "structCharInt" and "structCharChar"
        extractIndices <- function(type) {
            stdHeaders <- convention[convention$slots == slotName, ]$headers
            indices <- vapply(stdHeaders,
                              function(x) {
                                  if (structType == "structCharInt") {
                                      res <- unname(struct[names(struct) == x])
                                  } else if (structType == "structCharChar") {
                                      res <- which(dataHeaders == struct[names(struct) == x])
                                  } else {
                                      stop("Error with struct type.")
                                  }
                                  return(ifelse(length(res) == 0, 0, res))
                              }, numeric(1))
            return(indices)
        }
        # Main of the function as.indices
        if (structType == "structListInt") {
            return(struct[[slotName]])
        } else {
            return(extractIndices(structType))
        }
    }

    # Attribution ## Faire plutôt un : si non NULL alors attribution et data frame...
    idAll    <- seq_len(ncol(data))
    idSpace  <- as.indices("space", names(data), struct)
    idTime   <- as.indices("time", names(data), struct)
    idObs    <- as.indices("obs", names(data), struct)
    idLabels <- idAll[!(idAll %in% c(idSpace, idTime, idObs))] # Les NULL se retirent automatiquement avec c(...)

    labels <- as.data.frame(data[, idLabels])
    space  <- as.data.frame(data[, idSpace])
    time   <- as.data.frame(data[, idTime])
    obs    <- as.data.frame(data[, idObs])

    # Except for labels, force renaming to be sure the column names meet the convention
    names(labels) <- names(data)[idLabels]
    rename <- function(slotName, length) {
        return(convention[convention$slots == slotName, ]$headers[0:length])
    }
    names(space) <- rename("space", ncol(space))
    names(time)  <- rename("time", ncol(time))
    names(obs)   <- rename("obs", ncol(obs))

    if(any(dim(labels) == 0)) labels <- data.frame(NULL) # plutôt : si NULL.. alors
    if(any(dim(space) == 0))  space  <- data.frame(NULL) # plutôt : si NULL.. alors
    if(any(dim(time) == 0))   time   <- data.frame(NULL) # plutôt : si NULL.. alors
    if(any(dim(obs) == 0))    obs    <- data.frame(NULL) # plutôt : si NULL.. alors

    object <- structure(list(labels = labels,
                             space  = space,
                             time   = time,
                             obs    = obs),
                        class = c(type, "intensity", "list"))

    valid_intensity(object)

    object
}

#==============================================================================#
# Create intensity object
#==============================================================================#
#==============================================================================#
# Constructors
#==============================================================================#
#' count, incidence and severity objects
#'
#' Constructs an object of one of the three classes inheriting from
#' \code{intensity} class. The choice of the class depends on the nature of the
#' data sets.
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
#'
#' @param data A data frame containing all the data. Each line corresponding to a records.
#' @param struct A vector with all the corresponding variables. The different
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
#' @name intensity
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' @rdname intensity
#' @export
#------------------------------------------------------------------------------#
count <- function(data, mapping) { # Pas nécessairement besoin de x,y,z,t,r,n nécessaiere seuelment r et n ici pour incidence
    init_intensity(data, mapping, type = "count")
}

#------------------------------------------------------------------------------#
#' @rdname intensity
#' @export
#------------------------------------------------------------------------------#
incidence <- function(data, mapping) { # Pas nécessairement besoin de x,y,z,t,r,n nécessaiere seuelment r et n ici pour incidence
    init_intensity(data, mapping, type = "incidence")
}

#------------------------------------------------------------------------------#
#' @rdname intensity
#' @export
#------------------------------------------------------------------------------#
severity <- function(data, struct = NULL) { # Pas nécessairement besoin de x,y,z,t,r,n nécessaiere seuelment r et n ici pour incidence
    initIntensity(data, struct, type = "incidence")
}



#------------------------------------------------------------------------------#
#' @keywords internal
#------------------------------------------------------------------------------#
is.empty.field <- function(field) {
    stopifnot(is.data.frame(field))
    nrow(field) == 0
}


#==============================================================================#
# Subsetting
#==============================================================================#
# Extract parts of an \code{intensity} object
#
# Operators acting on \code{intensity} objects to extract parts.
#
# Conceptually, these operators work in the same way as the ones documented in
# the package \code{base}. Take a look at
# \code{base::}\code{\link[base]{Extract}} for further information.
#
# @param x Object from which to extract element(s).
# @param i,j,k Spatial indices (corresponding to \code{x}, \code{y} and
#   \code{z}, respectively) specifying elements to extract.
# @param t Time indice specifying elements to extract.
#
# @return \code{[} returns an \code{intensity} object, and \code{[[} returns a
# \code{data.frame}. \code{[[x]]} can be view as a shortcut for
# \code{as.data.frame(x)}.
#
# @name Extract
# @aliases [
#
# @export
#------------------------------------------------------------------------------#
#`[.intensity` <- function(x, i, j, k, t) {
#    data <- extractData(x, i, j, k, t)
#    initIntensity(data, type = class(x)[[1]]) # class(x)[[1]] => retourne la classe fille
#}

#------------------------------------------------------------------------------#
# @rdname Extract
# @name [[
# @export
#------------------------------------------------------------------------------#
#`[[.intensity` <- function(x, i, j, k, t) {
#    extractData(x, i, j, k, t)
#}

#------------------------------------------------------------------------------#
# @rdname Extract
# @name $
# @export
#------------------------------------------------------------------------------#
#`$.intensity` <- function(x, name) {
#    res <- which(str(x, simplify = TRUE, legacy = FALSE)$name == name)
#    if (length(res) == 0) stop("name does not exist.")
#    if (length(res) > 1)  stop("name is duplicated.") # Not supposed to occur
#    as.data.frame(x)[[name]]
#}

#------------------------------------------------------------------------------#
#' @keywords internal
#------------------------------------------------------------------------------#
######## To go ahead -----------------------------------------------------------------------
extractData <- function(x, i, j, k, t) {

    #ncol(x@space)
    if (missing(i))    i <- seq_len(nrow(x$space))
    if (missing(j))    j <- seq_len(nrow(x$space))
    #if (missing(k))    k <- # 3rd spatial dimension (i.e. z)... if any

    #ncol(x@time)
    if (missing(t))    t <- seq_len(nrow(x$time))

    idxSpace <- which((x$space$x %in% i) & (x$space$y %in% j)) # s'occuper du z !!!
    idxTime  <- which(x$time$t %in% t)
    idx      <- intersect(idxSpace, idxTime)
    if (length(idx) == 0) return(NA)

    space <- x$space[idx, , drop = FALSE]
    time  <- x$time[idx, , drop = FALSE]
    obs   <- x$obs[idx, , drop = FALSE]

    data.frame(space, time, obs)
}


#==============================================================================#
# is tools: intensity
#==============================================================================#
# Add a small doc here
#' @export
#------------------------------------------------------------------------------#
is.intensity <- function(object) return(is(object, "intensity"))
#' @export
is.count     <- function(object) return(is(object, "count"))
#' @export
is.incidence <- function(object) return(is(object, "incidence"))
#' @export
is.severity  <- function(object) return(is(object, "severity"))

#==============================================================================#
# as tools: incidence
#==============================================================================#
#' as.incidence & is.incidence
#'
#' Functions to check if an object is an \code{\link{incidence}} object, or coerce it if possible.
#'
#' \code{as.incidence} is a generic function with different methods. In this package, only
#' methods for \code{\link{count}} and \code{\link{severity}} objects are implemented.
#'
#' \code{is.incidence} is basically a wrapper.
#'
#' @param object Any R object.
#'
#' @return
#' \code{as.incidence} returns an \code{incidence} object if it is possible.
#'
#' \code{is.incidence} returns \code{TRUE} if its argument is an \code{incidence} object
#' (that is, has \code{incidence} amongst its classes) and \code{FALSE} otherwise.
#'
#' @seealso \code{\link{incidence}}
#' @name as.incidence
#------------------------------------------------------------------------------#
as.incidence <- function(object, ...) UseMethod("as.incidence")

#------------------------------------------------------------------------------#
#' @rdname as.incidence
#' @name as.incidence
#' @method as.incidence count
#' @export
#------------------------------------------------------------------------------#
as.incidence.count <- function(object, ...) {
    r <- as.numeric(object$obs$r > 0)
    n <- 1
    tmp <- data.frame(as.data.frame(object, fields = c("space", "time")),
                      r = r, n = n)
    count(tmp)
}

#------------------------------------------------------------------------------#
#' @rdname as.incidence
#' @name as.incidence
#' @method as.incidence severity
#' @export
#------------------------------------------------------------------------------#
as.incidence.severity <- function(object, ...) {
    stop("Not yet implemented.")
}


#==============================================================================#
# as tools: severity
#==============================================================================#

#------------------------------------------------------------------------------#
#' As severity
#'
#' TO DO
#'
#' @name as.severity
#' @export
#------------------------------------------------------------------------------#
as.severity <- function(object, ...) UseMethod("as.severity")

#------------------------------------------------------------------------------#
#' @method as.severity count
#' @export
#------------------------------------------------------------------------------#
as.severity.count <- function(object, ...) {
    stop("Not yet implemented.")
}

#------------------------------------------------------------------------------#
#' @method as.severity incidence
#' @export
#------------------------------------------------------------------------------#
as.severity.incidence <- function(object, ...) {
    r <- object$obs$r / object$obs$n
    tmp <- data.frame(as.data.frame(object, fields = c("space", "time")),
                      r = r)
    severity(tmp)
}


#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
str.incidence <- function(object, ..., legacy = TRUE, simplify = FALSE) {

    if (legacy) return(NextMethod())

    # Get structural information
    name    <- lapply(object, names)
    iMember <- lapply(name, function(x) seq_len(length(x)))
    i <- 0
    iAll <- lapply(iMember, function(x) {
        if (length(x) == 0) return(integer(0))
        indices <- (i + 1):(i + length(x))
        i <<- rev(indices)[1]
        indices
    })
    res <- Map(data.frame,
               name    = name,
               iMember = iMember,
               iAll    = iAll,
               stringsAsFactors = FALSE)

    if (simplify) {
        do.call(rbind, res)
    } else {
        res
    }
}


#==============================================================================#
# as.data.frame
#==============================================================================#
#' Coerce to a data frame
#'
#' Functions to coerce an \code{intensity} object to data.frame class.
#'
#' @param x An \code{intensity} object.
#' @inheritParams base::as.data.frame
#' @param members Slots of \code{x} to use to build the data frame.
#'
#' @return A data frame.
#'
#' @seealso \code{\link{count}} \code{\link{incidence}} \code{\link{severity}}
#'
#' @examples
#' myData <- incidence(Cochran1936)
#' head(as.data.frame(myData))
#' head(as.data.frame(myData, fields = c("space", "obs")))
#'
#' @name as.data.frame
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' @rdname as.data.frame
#' @name as.data.frame
#' @method as.data.frame intensity
#' @export
#------------------------------------------------------------------------------#
as.data.frame.intensity <- function(x, ..., members = c("labels", "space",
                                                        "time", "obs")) { # faire en sorte partial match avec fields // field ???
    stopifnot(all(members %in% names(x)))
    struct <- str(x, legacy = FALSE)
    res    <- x[members]
    res    <- res[!vapply(res, function(field) nrow(field) == 0, logical(1))] # As all the slots are data frames of same nrow or nroe == 0
    struct <- struct[names(res)]
    res    <- as.data.frame(res, ...)
    names(res) <- do.call(c, lapply(struct, function(x) x[["name"]]))
    attr(res, "struct") <- struct
    res
}

# Vérifier que tous les élément de la liste == même nature

#==============================================================================#
# as.array
#==============================================================================#
#' Test
#'
#' Test
#'
#' @param x An intensity object.
#' @param dim Name or ID of the dimensions.
#'
#' @examples
#'
#' myData <- incidence(Cochran1936)
#' as.array(myData, dim = c(24, 60, 2))
#'
#' @name as.array
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' @rdname as.array
#' @name as.array
#' @method as.array intensity
#' @export
#------------------------------------------------------------------------------#
as.array.intensity <- function(x, dim, ...) {
    if (missing(dim)) {
        dim <- c("x", "y", "t") # PRoblème si t == 2 apres un split... passer t à 1 autmomatiquement???
    } else {
        # 1. La dimension ne peut pas être contenu dans @obs
        # 2. Doit être contenu dans "convention"
    }
    obj <- as.data.frame(x)
    obj <- obj %>% dplyr::arrange_(.dots = rev(dim))
    maxis <- vapply(dim, function(i) length(unique(obj[, i])), numeric(1)) # integer(1)
    #maxis <- vapply(dim, function(i) max(obj[, i]), numeric(1))
    res <- array(obj$r, dim = maxis)
    if (class(x)[[1]] == "incidence") {
        res <- list(r = res, n = array(obj$n, dim = maxis))
    }
    return(res)
}

#==============================================================================#
# Misc
#==============================================================================#
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

#' @method as.character count
#' @export
as.character.count <- function(x, ...) "<count object>"

#' @method as.character incidence
#' @export
as.character.incidence <- function(x, ...) "<incidence object>"

#' @method as.character severity
#' @export
as.character.severity <- function(x, ...) "<severity object>"


#==============================================================================#
# print
#==============================================================================#
# intensity
#' @export
#------------------------------------------------------------------------------#
print.intensity <- function(x, ...) {
    mapped_data <- lapply(x$mapping, eval, envir = x$data, enclos = NULL)
    n_row       <- nrow(x$data)
    if ("t" %in% names(mapped_data)) {
        n_time <- length(unique(mapped_data$t))
    }
    # if (nTime > 2) {
    #     if (length(unique(diff(x$time$t)) == 1))
    #         nat <- "regular"
    #     else
    #         nat <- "irregular"
    # }

    # xLen <- length(levels(as.factor(x$space$x))) ## Devrait déjà être un facteur !!!!
    # yLen <- length(levels(as.factor(x$space$y))) ## Devrait déjà être un facteur !!!!
    ### Et pour z ????????????

    # nSU    <- nrow(unique( data.frame(x = x$space$x,
    #                                   y = x$space$y) )) # Trop specifique
    # if (nSU == 0) {
    #     georeferenced <- FALSE
    #     if (nTime == 1) {
    #         nSU <- nRow
    #     } else {
    #         nSU <- "" # to complete
    #     }
    # } else {
    #     georeferenced <- TRUE
    # }
    # anyNA     <- ifelse(any(is.na(x$obs)), TRUE, FALSE) # traiter du cas des NA dans le validity de incidence
    # complete  <- is.completeArray(x)
    # anyLabels <- nrow(x$labels)

    # word <- function(x) {if(x) "Yes" else "No"}

    # georeferenced <- word(georeferenced)
    # noNA          <- word(!anyNA)
    # complete      <- word(complete)
    # anyLabels     <- word(anyLabels)

    # More fancy display
    #chars <- nchar(c(nSU, nTime, georeferenced, noNA, complete))
    #ltab <- max(chars) - chars
    #ltab <- vapply(ltab, function(x) paste0(rep(" ", x + 2), collapse = ""), character(1))

    # Display prepartation
    # cat("<", class(object)[1], " object>\n",
    #     ltab[1], nSU, " sampling unit", ifelse(nSU > 1, "s", ""), " (", xLen, " × ", yLen, ")\n",
    #     ltab[2], sep = "")
    # if (nTime == 0) cat("1 (implicit)")
    # if (nTime <= 2) cat(nTime)
    # if (nTime > 2)  cat(nTime, " (", nat, " time-step)", sep ="")
    # cat(" snapshot", ifelse(nTime > 1, "s", "") ,"\n",
    #     ltab[3], georeferenced, " explicitly spatialized (", ncol(object@space), " dim(s))\n",
    #     ltab[4], noNA, " NA-free\n",
    #     ltab[5], complete, " a complete array\n",
    #     sep = "")

    nSU <- 2
    xLen <- 2
    yLen <- 2
    nat <- "TODO"
    anyLabels <- "TODO"
    georeferenced <- "TODO"
    noNA <- complete <- "TODO"


    cat("<", class(x)[[1]], " object>\n",
        "\u2022 number of sampling unit", ifelse(nSU > 1, "s: ", ": "),
        nSU, " (", xLen, " × ", yLen, ")\n",
        "\u2022 number of snapshot", ifelse(n_time > 1, "s: ", ": "),
        sep = "")
    if      (n_time == 0) cat("1 (implicit)\n")
    else if (n_time <= 2) cat(n_time, "\n", sep = "")
    else                  cat(n_time, " (", nat, " time-step)\n", sep ="")
    cat("\u2022 any labels: ", anyLabels, "\n", sep = "")
    cat("\u2022 explicitly spatialized: ", georeferenced,
        " (", "ncol(x$space)", " dim(s))\n",
        "\u2022 NA-free: ", noNA, "\n",
        "\u2022 complete array: ", complete, "\n", sep = "")

    ## EXPORTER TOUTES LES METHODS QUI PERMETTENT l'AFFICHAGE DE CES INDICATEURS
    ## Toolbox:
    # nSU    # Number of sampling units (if same coordinates over time = same sampling units)
    # Dire si homogene au cours temps, ou nom
    # Le mieux, retourner un S3
    # ... <- nTime() # Besoin de savoir
    # structure(list(), class = "NSU")
    # nSpace #
    # nTime  # number of temporal snapshots
    # is.explicitSpatialized #
    # is.NAfree              # ~~~
    # is.completeArray       #
}


# <incidence object>
# - number of sampling units : 1440 (24 × 60) -> nSU() => return: list(total = ..., detail = list(x = , y = ,...))
# - number of snapshots -----: 2 -------------------> nTimes() => return: list(total = ...)
# - any labels --------------: No ---------------------------> is.anyLabel()
# - explicitly spatialized --: Yes (2 dim(s)) ---> is.explicitSpace()
# - NA-free -----------------: Yes -----------------------------> is.anyNA()
# - complete array ----------: Yes ----------------------> is.completeArray(): OK


## summary(object, ..., rescale = FALSE)
#
# General statistics:
#
# time median mean   sd   se
#    1    4.3  4.5  ...  ...
#    2    5.2  5.2  ...  ...
#
# Indexes of aggregation:
#
# time: 1
# Name           Value  test  stat   p-value
# - Fisher's     1.52   C(a)  23.07  0.06    *
# - Morisita's   1.21   .     .      .
#
# time: 2
# Name           Value  test  stat   p-value
# - Fisher's     1.52   C(a)  23.07  0.06    *
# - Morisita's   1.21   .     .      .


#' @export
is.completeArray <- function(object) {
    dt <- as.data.frame(object, fields = c("space", "time"))
    ref <- expand.grid(lapply(dt, unique))
    #names(dt) <- names(ref) <- make.names(rep("X", ncol(dt)), unique = TRUE)
    res <- merge(dt, ref, all = TRUE)
    if (nrow(res) == nrow(dt)) return(TRUE)
    else return(FALSE)
}


#==============================================================================#
# plot
#==============================================================================#

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
plot.count <- function(x, y, ...) {
    # if time ???
    mapped_data <- map_data(x)
    with(mapped_data, plot(x, y, bg = gray(1 - (r / max(r))), cex = 4, pch = 21))
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
plot.intensity <- function(x, y, type = c("map_2D", "progress_curve"), t = NULL,
                           nRibbon = 2, ...) {
    g <- list()
    if (!is.null(t)) {
        if ("map_2D" %in% type)
            g[[length(g) + 1]] <- plotmap_2D(x, t, nRibbon)
        if ("progress_curve" %in% type)
            g[[length(g) + 1]] <- plotprogress_curve(x, t)
        ### to correct when t only 1 value
        #geom_path: Each group consist of only one observation. Do you need to adjust the group aesthetic?
    } else {
        timeSteps <- unique(x$time$t)
        if ("map_2D" %in% type)
            g[[length(g) + 1]] <- plotmap_2D(x, timeSteps, nRibbon)
        if ("progress_curve" %in% type)
            g[[length(g) + 1]] <- plotprogress_curve(x, timeSteps) # to generaliser
    }
    for (i in seq_len(length(g)))
        print(ggplot() + g[[i]])
}

#------------------------------------------------------------------------------#
#' @keywords internal
#------------------------------------------------------------------------------#
plotmap_2D <- function(object, time, nRibbon = 2) { # Proposer l'inversion visuelle x <-> y
    if (!is(object, "intensity")) stop("object must be an intensity object.")

    data <- subset(as.data.frame(object), t %in% time)

    xScale   <- 1:max(data$x) # pas top ; tmp
    yScale   <- 1:max(data$y) # pas top ; tmp
    if (class(object)[[1]] == "incidence")
        maxScale <- max(data$n) # Gérer le cas différents n dans même jeu de données ; marche seulement pour incidence !!!
    else
        maxScale <- max(data$r)

    #library("ggplot2")
    #Cochran1936$texp <- with(Cochran1936, factor(t, labels = c("18 December 1929",
    #                                                           "31 December 1929")))
    #g <- ggplot(Cochran1936, aes(x = x, y = y, colour = as.factor(r)))
    #g <- g + geom_point()
    #g <- g + scale_colour_manual(values = c("white", "red"))
    #g <- g + facet_grid(. ~ texp) + theme(legend.position="none")
    #g <- g + coord_flip() + scale_x_reverse()
    #print(g)

    g <- list(
        #----geom_point(data = data, aes(x = x, y = y, color = r), size = 10, shape = 15),
        geom_tile(data = data, aes(x = x, y = y, fill = r)),
        ## Valeurs fixées à généraliser ci-après !!!!!!!!
        #scale_colour_gradient(name = "Disease\nintensity",
        #                      low = "white", high = "red", breaks = seq(0, maxScale),
        #                      limits = c(0, maxScale), guide = "legend"),
        scale_fill_gradient(low = "white", high = "red",
                            breaks = seq(0, maxScale),
                            limits = c(0, maxScale), guide = "legend"),
        # Faire un conditionel ici, de la forme:
        # si pas "reverse x", alors:

        scale_x_continuous(breaks = xScale, expand = c(0,0)),
        scale_y_continuous(breaks = yScale, expand = c(0,0)),

        # sinon
        #scale_x_reverse(breaks = xScale, expand = c(0,0)),
        #

        #scale_y_discrete(breaks = yScale, expand = c(0,0)),

        #                         expand = c(0,0)),
        ##
        ## expand = c(0,0) : - OK avec geom_tile
        ##                   - NON avec geom_point
        ##

        theme(panel.grid = element_blank()),

        #coord_flip(), # Le faire en conditionnel, avec formula : x ~ y (default), ou y ~ x
        coord_fixed(), # Pour avoir des carr'es pour sûr
        facet_wrap(~ t, nrow = nRibbon)#, # to wrap a 1d ribbon of panels into 2d.
        ### to generaliser pour nrow et ncol

        ### [BEG] TEMPORARY - FOR 2015 APS ANNUAL MEETING
        #labs(x = "", y = ""),
        #scale_x_reverse(breaks = NULL, labels = NULL),
        #scale_y_discrete(breaks = NULL, labels = NULL)
        ### [END] TEMPORARY
    )
    #g <- g + labs(title = paste0("Field: ",nom.feuille,"\n","Date: ",sub.sub.data$date[1]))
    #ggsave !!!! c'est une bonne idée de proposer le choix !!!
    ### + proposer l'extrapolation :
    # library(akima)
    # x <- seq(1,10, length=200)
    # y <- seq(1,20, length=200)
    # data2 <- interp(sub.sub.data$row, sub.sub.data$sample, sub.sub.data$dleaflet, xo=x, yo=y, linear=TRUE, extrap=FALSE)
    # data2 <- data.frame(expand.grid(row=data2$x, sample=data2$y), dleaflet=melt(data2$z)$value)
    #
    return(g)
}

#------------------------------------------------------------------------------#
#' @keywords internal
#------------------------------------------------------------------------------#
plotprogress_curve <- function(object, time, col = "red", visuRecap = "all") {
    if (!is(object, "intensity")) stop("object must be an intensity object.")
    #i <- max(sub.data$scoring)
    data  <- subset(as.data.frame(object), t %in% time)
    recap <- data %>%
        dplyr::group_by(t) %>%
        dplyr::summarise(
            meanD = mean(r),
            sdD = sd(r))
    #g <- ggplot(data = recap, aes(x = t, y = meanD))
    #g <- ggplot()
    g <- list(
        scale_y_continuous("score (1-9 scale)", breaks=0:9), #### à Généralsier !!!!
        ######### De plus attention, warning quand on combine les graphs
        #### Scale for 'y' is already present. Adding another scale for 'y', which will replace the existing scale.
        #g <- g + scale_x_continuous("date", breaks=seq(1,i))
        geom_point(data = data, aes(x = t, y = r),
                   alpha = 0.2, size = 3, colour = col,
                   position = position_jitter(w = 0.2, h = 0.2)),
        geom_line(data = recap, aes(x = t, y = meanD),
                  size = 1.5, linetype = 2, colour = col),
        geom_point(data = recap, aes(x = t, y = meanD),
                   size = 6, shape = 15, colour = col)
    )

    #g <- g + labs(title = paste0("ALS score over time","\n","Field: ",nom.feuille))
    #ggsave(filename=paste0("fig/",nom.feuille,"-over-time.png"), plot = g)
    if (visuRecap == "all")
        g <- g[c(TRUE, TRUE, TRUE, TRUE)] # Faire un truc plus intelligent ----------
    if (visuRecap == "line")
        g <- g[c(TRUE, TRUE, TRUE, FALSE)] # Faire un truc plus intelligent ----------
    if (visuRecap == "point")
        g <- g[c(TRUE, TRUE, FALSE, TRUE)] # Faire un truc plus intelligent ----------

    return(g)
}


#==============================================================================#
# get and set internal strure
#==============================================================================#
# Useful
#
# Useful
#
# @param object An object inherited from class \code{intensity}.
# @param simplify Is the return value stacked as an unique data frame?
# @param value Must be of the folowing form. Not possible for \code{obs} slot.
# @param legacy Legacy behavior of the \code{\link[utils]{str}} function.
#
# @examples
# \dontrun{
# myData <- incidence(Cochran1936)
# str(myData)
# str(myData) <- list(label = ,
#                        space = c(x = 2, y = 1),
#                        time = c(t = 3))
#
# # To reverse x and y values:
# myData <- incidence(Cochran1936)
# plot(myData, type = "map_2D")
# str(myData) <- c(x = 2, y = 1)
# plot(myData, type = "map_2D")
# }
#
# @name str
#------------------------------------------------------------------------------#
#NULL

#------------------------------------------------------------------------------#
# @name str
# @rdname str
# @export
#------------------------------------------------------------------------------#
# str.intensity <- function(object, ..., simplify = TRUE, legacy = FALSE) {
#     if (legacy) {
#         return(utils::str(object, ...))
#     }
#     else {
#         i <- 0
#         slotLabels <- slotNames(object)
#         res <- lapply(seq_len(length(slotLabels)), function(j) {
#             data     <- slot(object, slotLabels[j])
#             if (ncol(data) == 0) return(NA)
#             colNames <- colnames(data)
#             indices  <- (i + 1):(i + length(colNames))
#             i <<- i + length(colNames)
#             return(data.frame(name = colNames,
#                               index = indices,
#                               summary = I(lapply(data, summary)),
#                               stringsAsFactors = FALSE))
#         })
#         names(res) <- slotLabels
#         if (simplify)
#             res <- do.call(rbind.data.frame, lapply(res, identity))
#         return(res)
#     }
# }


#------------------------------------------------------------------------------#
# @rdname str
# @name str<-
# @export
#------------------------------------------------------------------------------#
#setGeneric("str<-", function(object, value){
#    standardGeneric("str<-")
#})

#------------------------------------------------------------------------------#
# @rdname str
# @name str<-
# @export
#------------------------------------------------------------------------------#
#setMethod("str<-",
#          signature(object = "intensity"),
#          function(object, value) {
#              # Uniquement pour str(object, simplify = TRUE) maintenant
#              # Pas sur que le passage par un as.data.frame soit vreiment efficace et sur (pour les labels)
#              ref <- str(object, simplify = TRUE)
#              tmp <- as.data.frame(object)
#              lapply(seq_len(length(value)), function(i) {
#                  names(tmp)[value[i]] <<- names(value)[i]
#              })
#              return(new(class(object), tmp))
#          }
#)

### S'occuper des labels correctement !!!!!!


#==============================================================================#
# Split and unsplit
#==============================================================================#
#' Divide into groups
#'
#' Great function
#'
#' @param x intensity object containing values to be divided into groups.
#' @param f A \code{factor} in the sense that \code{\link{as.factor}(f)} defines the grouping.
#'
#' @name split
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' @rdname split
#'
#' @examples
#' myData <- incidence(Cochran1936)
#'
#' res1 <- split(myData, myData@time[, 1])
#' # is equivalent to:
#' res2 <- lapply(levels(as.factor(myData@time[, 1])),
#'                function(i) myData[t = i])
#' names(res2) <- levels(as.factor(myData@time[, 1]))
#'
#' @export
#------------------------------------------------------------------------------#
split.incidence <- function(x, f, drop, ...) { # S'occuper du drop !!!
    nameVars <- lapply(x, names) # À utiliser après
    data <- as.data.frame(x)
    data <- split(data, f, drop, ...)
    lapply(data, function(sub) initIntensity(sub, type = class(x)[[1]])) # Que se passe-t-il avec les labels ?, on utilise pas struct ici???
}

#------------------------------------------------------------------------------#
#' Super fonction
#'
#' @name unsplit
#' @export
#------------------------------------------------------------------------------#
unsplit <- function(value, f, drop = FALSE) UseMethod("unsplit")

#------------------------------------------------------------------------------#
#' @rdname unsplit
#' @export
#------------------------------------------------------------------------------#
unsplit.list <- function(value, f, drop) { # S'occuper du drop !!! Se transmet-il de la g'en'eric `a ici?

    if(!all(vapply(value, is.intensity, logical(1)))) # on envoie vers le unsplit standart !!
        unsplit(value, f, drop)
    ## +++ Plus tous identique !!! m^eme sous-classe
    ## Headers identiques !!!!
    ##### NOTE : f not used hereinafter
    class    <- class(value[[1L]])
    #headers  <- lapply(str(value[[1L]], simplify = FALSE), function(i) i$index)
    headers  <- lapply(str(value[[1L]], simplify = FALSE), function(i) {
        if (!any(is.na(i))) i$index else NULL})
    headers <- headers[!sapply(headers, is.null)]

    value    <- lapply(value, as.data.frame)
    value    <- do.call(rbind.data.frame, ### Attention aux factor !!!!!
                        unname(value))    ### ici on perds de l'information avec unname... pas bien
    #df  <- lapply(value, as.data.frame)
    #res <- unsplit(value, f, drop) # S'occuper du drop !!! Se transmet-il de la g'en'eric `a ici?

    initIntensity(value, struct = headers, type = class(x)[[1]])

    #df <- do.call(rbind.data.frame, lapply(seq_len(length(value)), function(i) {
    #    res <- as.data.frame(value[[i]])
    #    res <- cbind.data.frame(f = f[i], res) # continue de faire des factors en th'eorie (v'erifier comment contre-carrer cela)
    #    return(res)}))
    #return(new(class(value[[1]]), df, struct = c(
    #    x = "x", y = "y", z = "z", t = "f", r = "r", n = "n"
    #)))#, struct = headers)) ### SUPER TMP !!!!!!!!!!!!!!!!!!!!!11
}

#------------------------------------------------------------------------------#
#' regroup data
#'
#' Useful.
#'
#' @param unitDim Dimensions of the base unit. All the dimensions must be an
#' integer diviser of the size of the object. Basically should be a list of two named
#' element: space and time, each one speciffying the dimensions. If not a list,
#' must be the same size as \code{length(object@@space, object@@time)}.
#' @param groupBy Parameter in the form of c("", "", ...). A vector of names
#' @param fun The function to use.
#'
#' @examples
#' myData <- incidence(Cochran1936)
#' myData <- regroup(myData, unitDim = c(3, 3, 1))
#' myData <- as.severity(myData)
#'
#' ## The folowing code do not work because it leads to non-integer number
#' ## which is imcomatible with the definition of incidence
#' \dontrun{
#' myData <- regroup(myData, unitDim = c(3, 3, 1), fun = mean)
#' }
#' ## To perform the desired calculation, we need to do:
#' ## Because we go out of the strict definition of the class incidence
#' myData <- incidence(Cochran1936)
#' myData <- regroup(myData, unitDim = c(3, 3, 1))
#' ## With dplyr library:
#' as.data.frame(myData) %>% dplyr::summarise(mean = r / n)
#' ## Or with a vanilla synthax:
#' myData <- as.data.frame(myData)
#' myData$mean <- myData$r / myData$n
#' @export
#------------------------------------------------------------------------------#
regroup <- function(object, ...) UseMethod("regroup")

#------------------------------------------------------------------------------#
#' @rdname regroup
#' @export
#------------------------------------------------------------------------------#
regroup.intensity <- function(object, unitDim, groupBy, fun = sum,
                              keepLabels = TRUE) {

    # Pour unitDim attention, doit respecter certaines contraines:
    # competer grid, les nombres se suivent, des entier de pr'ef'erence..
    # beaucoup de limites en fait !!!

    if (!missing(unitDim) & !missing(groupBy))
        stop("unitDim and groupBy cannot be specified together.")
    if (!missing(unitDim)) {

        ## Gérer le keepLabels ICI !!!

        # Ci-après, il manque bp de vérifications !!!!
        if (length(unitDim) > 2)
            stop("Only size <= 2 taken into account so far.")
        if (any(unitDim < 1))
            stop("unitDim sizes must be >= 1.")
        #if (modulo .... )
        # stop(pas un multiple entier !!!!!)
        # ci-après ne marche que si x et y sont des id de cases (entier, case sans trous)
        #space[,1] if any !!!!!
        object$space[, 1] <- ceiling(object$space[, 1] / unitDim[1])
        #space[,2] if any !!!!!
        object$space[, 2] <- ceiling(object$space[, 2] / unitDim[2])

        obsNames <- str(object)$obs$name

        dots <- list(~fun(r), ~fun(n))

        newObject <- as.data.frame(object) %>%
            dplyr::group_by(x, y, t) %>% ## Tous les noms des labels + space + time
            #dplyr::group_by_(names/str(object@labels), x, y, t) %>% ## Tous les noms des labels + space + time
            # Attention: " Error: unknown variable to group by : t" si t n'existe pas !!!!
            # "epiphy::regroup(unit = c(y=1, x=10))" marche pas snif.
            dplyr::summarise_(.dots = setNames(dots, obsNames))

        headers  <- lapply(str(object, simplify = FALSE), function(i) {
            if (!any(is.na(i))) i$index else NULL})
        headers <- headers[!sapply(headers, is.null)]
        #browser()
        return(new(class(object),
                   newObject,
                   struct = headers))
    } else if (!missing(groupBy)) {
        if (keepLabels)
            dots1 <- lapply(c(str(object, simplify = FALSE)$label$name, groupBy), as.symbol)
        else
            dots1 <- lapply(groupBy, as.symbol)

        obsNames <- str(object, simplify = FALSE)$obs$name
        dots <- list(~fun(r), ~fun(n))

        newObject <- as.data.frame(object) %>%
            dplyr::group_by_(.dots = dots1) %>%
            dplyr::summarise_(.dots = setNames(dots, obsNames))

        #headers2  <- lapply(str(object, simplify = FALSE), function(i) i$index)

        return(initIntensity(newObject, type = class(x)[[1]]))
        #                             struct = headers2))
    } else {
        stop("unitDim or groupBy must be specified.")
    }

}



# @export
#dplyr::`%>%`

#------------------------------------------------------------------------------#
# @export
#------------------------------------------------------------------------------#
#setGeneric("%df>%", function(lhs, rhs) {
#    standardGeneric("%df>%")
#})

#------------------------------------------------------------------------------#
# @export
#------------------------------------------------------------------------------#
#setMethod("%df>%",
#       signature(lhs = "intensity"),
#       function(lhs, rhs) {
#           call <- match.call()
#           call$lhs <- quote(as.data.frame(lhs))
#           call[[1L]] <- as.name("%>%")
#           eval(call) # lhs est evalué dans cet environnement
#       }
#)



#member
#name
#index

