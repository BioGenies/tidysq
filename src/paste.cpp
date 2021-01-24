#include "tidysq.h"
#include "tidysq/paste.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_paste(const Rcpp::List& list_of_x,
                     const tidysq::Letter &NA_letter) {
    std::vector<Sq<RCPP_IT>> list_of_sq;
    for (const Rcpp::internal::const_generic_proxy<19, Rcpp::PreserveStorage>& x : list_of_x) {
        list_of_sq.push_back(import_sq_from_R(x, NA_letter));
    }
    return export_to_R(paste<RCPP_IT>(list_of_sq));
}
