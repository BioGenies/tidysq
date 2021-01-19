#include <Rcpp.h>

#include "tidysq/ops/complement.h"
#include "tidysq/Rcpp-import.h"
#include "tidysq/Rcpp-export.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_complement(const Rcpp::List& x,
                          const tidysq::Letter &NA_letter) {
    return export_to_R(complement<RCPP_IT>(import_sq_from_R(x, NA_letter)));
}
