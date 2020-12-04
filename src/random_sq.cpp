#include <Rcpp.h>

#include "tidysq/ops/random_sq.h"
#include "tidysq/Rcpp-import.h"
#include "tidysq/Rcpp-export.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_random_sq(const Rcpp::IntegerVector &n,
                         const Rcpp::IntegerVector &len,
                         const Rcpp::StringVector &alphabet,
                         const Rcpp::LogicalVector &use_gap) {
//    srand(seed);
    if (len.length() == 1) {
        return export_to_R(random_sq<RCPP_IT>(
                util::convert_to_scalar(n),
                util::convert_to_scalar(len),
                import_alphabet_from_R(alphabet, constants::DEFAULT_NA_LETTER),
                util::convert_to_scalar(use_gap)));
    }
    return export_to_R(random_sq<RCPP_IT>(
            Rcpp::as<std::vector<LenSq>>(len),
            import_alphabet_from_R(alphabet, constants::DEFAULT_NA_LETTER),
            util::convert_to_scalar(use_gap)));
}
