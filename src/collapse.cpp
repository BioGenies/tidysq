#include "tidysq/Rcpp-import.h"
#include "tidysq/Rcpp-export.h"
#include "tidysq/collapse.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_collapse(const Rcpp::List& x,
                        const tidysq::Letter &NA_letter) {
    return export_to_R(collapse<RCPP_IT>(import_sq_from_R(x, NA_letter)));
}
