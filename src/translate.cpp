#include "tidysq/Rcpp-import.h"
#include "tidysq/Rcpp-export.h"
#include "tidysq/ops/translate.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_translate(const Rcpp::List &x,
                         const int &table,
                         const tidysq::Letter &NA_letter) {
    return export_to_R(translate<RCPP_IT>(
            import_sq_from_R(x, NA_letter),
            table));
}
