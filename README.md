[![Build Status](https://api.travis-ci.org/michbur/tidysq.png)](https://travis-ci.org/michbur/tidysq)
[![codecov.io](https://codecov.io/github/michbur/tidysq/coverage.svg?branch=master)](https://codecov.io/github/michbur/tidysq?branch=master) 


## tidysq package

This package contains tools for analysis of biological sequences of amino acids (e.g., peptides, proteins) or nucleic acids (e.g., RNA, DNA). Through the efficient compression of the sequence data, *tidysq* allows studies of very large datasets in **R**.

### Installation

You can install the latest development version of the package using the `devtools` R package.

```
source("https://install-github.me/michbur/tidysq")
```

## Troubleshooting

This package heavily utilizes Rcpp to assure the fastest compression of the sequence data. If you have any compiler-related issues, as "clang: error: unsupported option '-fopenmp'", please follow [these instructions](https://github.com/RcppCore/RcppArmadillo/issues/143).

## Citation

For citation type:

```
citation("tidysq")
```

or use:
Michal Burdukiewicz, Dominik Rafacz, Weronika Puchala, Filip Pietluch, Katarzyna Sidorczuk, Stefan Roediger and Leon Eyrich Jessen (2019). tidysq: N-Gram Analysis of Biological Sequences. R package version 1.0. 
