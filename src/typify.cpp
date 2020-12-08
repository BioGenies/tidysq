#include <Rcpp.h>

#include "tidysq/Rcpp-import.h"
#include "tidysq/Rcpp-export.h"
#include "tidysq/ops/typify.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_typify(const Rcpp::List& x,
                      const Rcpp::StringVector& dest_type,
                      const Rcpp::StringVector& NA_letter) {
    return export_to_R(typify<RCPP_IT>(import_from_R(x, NA_letter),
            util::sq_type_for_sq_type_abbr(util::convert_to_scalar(dest_type))));
}
