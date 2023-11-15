## Resubmission

This is a resubmission. In this version I have:

* added a \value section for all exported functions where such section was missing.
    * clump.Rd
    * split.intensity.Rd
    * threshold.Rd

The package was archived on CRAN because of the use of 'default.stringsAsFactors'. See:
https://cran-archive.r-project.org/web/checks/2022/2022-03-10_check_results_epiphy.html
Uses of 'default.stringsAsFactors' have been removed to make the package compatible with R>=4.

## Test environments

* local Linux Mint 21.2, R 4.3.2
* GitHub Actions:
    * macos-latest (release), R 4.3.2
    * windows-latest (release), R 4.3.2
    * ubuntu-latest (devel), R-devel
    * ubunut-latest (release), R 4.3.2
    * ubunut-latest (oldrel-1), R 4.2.3
* win-builder:
    * release, R 4.3.2
    * devel, R-devel
* macOS builder, R 4.3.0

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 1 note ✖

* There was only 1 NOTE (this is a new release):
    * Maintainer: 'Christophe Gigot <ch.gigot@gmail.com>'
    * New submission
    * Package was archived on CRAN

## Reverse dependencies

This are no reverse dependencies.

