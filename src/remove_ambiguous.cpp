#include <Rcpp.h>

#include "tidysq.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_remove_ambiguous(const Rcpp::List& x,
                                const bool &by_letter,
                                const tidysq::Letter& NA_letter) {
    return export_to_R(remove_ambiguous<RCPP_IT>(import_from_R(x, NA_letter), by_letter));
}
