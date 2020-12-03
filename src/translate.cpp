#include "tidysq/Rcpp-import.h"
#include "tidysq/Rcpp-export.h"
#include "tidysq/ops/translate.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_translate(const Rcpp::List &x,
                         const Rcpp::IntegerVector &table,
                         const Rcpp::StringVector &NA_letter,
                         const Rcpp::LogicalVector &interpret_as_stop) {
    // TODO: replace with Rcpp::IntegerVector and coercing to scalar int
    return export_to_R(translate<RCPP_IT>(
            import_from_R(x, NA_letter),
            util::convert_to_scalar(table),
            util::convert_to_scalar(interpret_as_stop)));
}
