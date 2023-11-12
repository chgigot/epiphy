## Submission

In this submission, we have:

* removed uses of function 'default.stringsAsFactors' to make the package compatible with R>=4.
* added a 'return value' section for all exported functions where such section was missing.
* added 'cph' and 'ctb' roles in DESCRIPTION.

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

