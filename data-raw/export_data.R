# To test the SADIE procedure:

#------------------------------------------------------------------------------#
# Export aphid counts (Perry et al., 1999)
#------------------------------------------------------------------------------#
aphid_counts <- read.csv("data-raw/aphid_counts.csv", header = TRUE)

with(aphid_counts, plot(x, y, pch = NA))
with(aphid_counts, text(x, y, labels = d))

devtools::use_data(aphid_counts)

#------------------------------------------------------------------------------#
# Export arthropods counts (Holland et al., 1999)
#------------------------------------------------------------------------------#
arthropods_counts <- read.csv("data-raw/arthropods_counts.csv")
arthropods_counts_splitted <- split(arthropods_counts, arthropods_counts$set)

opar <- par()
par(mar = c(2,2,2,2)) # We do not see x and y labels anymore.
layout(matrix(1:6, nrow = 3, ncol = 2, byrow = TRUE))
invisible(lapply(arthropods_counts_splitted, function(set) {
    with(set, plot(x, y, pch = NA))
    with(set, text(x, y, labels = d))
}))
par(opar)
layout(1)

devtools::use_data(arthropods_counts)

#------------------------------------------------------------------------------#
# Export codling moth counts (Lavigne et al., 2010)
#------------------------------------------------------------------------------#

