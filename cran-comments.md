## Test environments
* local R installation, R 4.0.4
* ubuntu 16.04 (on travis-ci), R 4.0.4
* win-builder (devel)
https://github.com/BioGenies/tidysq/actions

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

The package contains a large amount of Cpp code which extends the building time. We weren't able to find any guidelines regarding the compilation time so we hope it is fine to submit the package as is.