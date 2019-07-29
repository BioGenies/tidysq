#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::IntegerVector integerize(Rcpp::RawVector raw_vec) {
  Rcpp::IntegerVector ret(raw_vec);
  return ret;
}

// [[Rcpp::export]]
Rcpp::RawVector rawize(Rcpp::IntegerVector int_vec) {
  Rcpp::RawVector ret(int_vec);
  return ret;
}
