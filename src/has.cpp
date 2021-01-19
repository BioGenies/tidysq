#include "tidysq/Rcpp-import.h"
#include "tidysq/Rcpp-export.h"
#include "tidysq/ops/has.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::LogicalVector CPP_has(const Rcpp::List &x,
                            const std::vector<std::string> &motifs,
                            const tidysq::Letter &NA_letter) {
   return Rcpp::wrap(has<RCPP_IT>(import_sq_from_R(x, NA_letter), motifs));
}
