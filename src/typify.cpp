#include <Rcpp.h>

#include "tidysq.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_typify(const Rcpp::List& x,
                      const std::string& dest_type,
                      const tidysq::Letter& NA_letter) {
    return export_to_R(typify<RCPP_IT>(import_from_R(x, NA_letter),
            util::sq_type_for_sq_type_abbr(dest_type)));
}
