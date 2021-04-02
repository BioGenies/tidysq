## Test environments
* local R installation, R 4.0.4
* ubuntu 16.04 (on travis-ci), R 4.0.4
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

Updated tests to be compliant with testthat release 3. 

Removed explicit declarations of default move constructors, copy constructors, move assignment operators and copy assignment operators which may cause problems on some platforms.