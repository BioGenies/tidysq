#include <Rcpp.h>

#include "tidysq/exports.h"
#include "tidysq/remove_ambiguous.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_remove_ambiguous(const Rcpp::List& x,
                                const bool by_letter,
                                const Rcpp::StringVector& NA_letter) {
    return export_to_R(remove_ambiguous<RCPP>(import_from_R(x, NA_letter), by_letter));
}
