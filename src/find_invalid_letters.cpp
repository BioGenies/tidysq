#include "tidysq/ops/find_invalid_letters.h"
#include "tidysq/Rcpp-import.h"
#include "tidysq/Rcpp-export.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_find_invalid_letters(const Rcpp::List& x,
                                    const std::string& dest_type,
                                    const tidysq::Letter& NA_letter) {
    return Rcpp::wrap(find_invalid_letters(import_sq_from_R(x, NA_letter),
                                           util::sq_type_for_sq_type_abbr(dest_type)));
}
