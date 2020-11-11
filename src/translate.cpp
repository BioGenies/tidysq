#include <Rcpp.h>

#include "tidysq/exports.h"
#include "tidysq/translate.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_translate(const Rcpp::List &x,
                         const int &table,
                         const Rcpp::StringVector &NA_letter) {
    // TODO: replace with Rcpp::IntegerVector and coercing to scalar int
    return export_to_R(translate<RCPP>(import_from_R(x, NA_letter), table));
}
