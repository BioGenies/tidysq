#include "tidysq/ops/bite.h"
#include "tidysq/ops/skip.h"
#include "tidysq/Rcpp-import.h"
#include "tidysq/Rcpp-export.h"
#include "tidysq/util/handle_warning_message.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_bite(const Rcpp::List& x,
                    const Rcpp::IntegerVector& indices,
                    const tidysq::Letter &NA_letter,
                    const std::string &on_warning) {
    if (Rcpp::is_true(Rcpp::all(indices > 0))) {
        Rcpp::IntegerVector cpp_indices = indices - 1;
        util::ResultWrapper<Sq<RCPP_IT>> ret =
                bite<RCPP_IT>(import_sq_from_R(x, NA_letter), Rcpp::as<std::vector<LenSq>>(cpp_indices));
        if (ret.has_message()) {
            util::handle_warning_message(ret.message_text(), on_warning);
        }
        return export_to_R(ret.result());
    } else if (Rcpp::is_true(Rcpp::all(indices < 0))) {
        Rcpp::IntegerVector cpp_indices = -indices - 1;
        Sq<RCPP_IT> ret = skip<RCPP_IT>(import_sq_from_R(x, NA_letter), Rcpp::as<std::vector<LenSq>>(cpp_indices));
        return export_to_R(ret);
    } else {
        throw std::invalid_argument("indices must be either all positive or all negative");
    }
}
