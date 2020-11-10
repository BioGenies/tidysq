#include <Rcpp.h>
#include <iterator>

#include "tidysq/exports.h"
#include "tidysq/find_invalid_letters.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_find_invalid_letters(const Rcpp::List& x,
                                    const Rcpp::StringVector& dest_type,
                                    const Rcpp::StringVector& NA_letter) {
    return Rcpp::wrap(find_invalid_letters(import_from_R(x, NA_letter),
            util::sq_type_for_sq_type_abbr(util::get_scalar_string_value(dest_type))));
}
