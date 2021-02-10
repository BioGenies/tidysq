#include <Rcpp.h>

#include "tidysq/ops/random_sq.h"
#include "tidysq/Rcpp-import.h"
#include "tidysq/Rcpp-export.h"

using namespace tidysq;

// [[Rcpp::export]]
Rcpp::List CPP_random_sq(const int &n,
                         const Rcpp::IntegerVector &len,
                         const Rcpp::StringVector &alphabet,
                         const bool &use_gap) {
//    srand(seed);
    if (len.size() == 1) {
        return export_to_R(random_sq<RCPP_IT>(
                n,
                util::convert_to_scalar(len),
                import_alphabet_from_R(alphabet),
                use_gap));
    }
    return export_to_R(random_sq<RCPP_IT>(
            Rcpp::as<std::vector<LenSq>>(len),
            import_alphabet_from_R(alphabet),
            use_gap));
}
