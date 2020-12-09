#include "tidysq/Rcpp-import.h"
#include "tidysq/Rcpp-export.h"
#include "tidysq/ops/has.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::LogicalVector CPP_has(const Rcpp::List &x,
                            const Rcpp::StringVector &motifs,
                            const Rcpp::StringVector &NA_letter) {
   return Rcpp::wrap(has<RCPP_IT>(import_from_R(x, NA_letter), util::convert_string_vector(motifs)));
}
