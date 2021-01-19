#include <Rcpp.h>

#include "tidysq/Rcpp-import.h"
#include "tidysq/Rcpp-export.h"
#include "tidysq/ops/remove_NA.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_remove_NA(const Rcpp::List& x,
                         const bool &by_letter,
                         const tidysq::Letter& NA_letter) {
    return export_to_R(remove_NA<RCPP_IT>(import_sq_from_R(x, NA_letter), by_letter));
}
