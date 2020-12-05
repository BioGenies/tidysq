#include <Rcpp.h>

#include "tidysq.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_translate(const Rcpp::List &x,
                         const int &table,
                         const tidysq::Letter &NA_letter,
                         const bool &interpret_as_stop) {
    return export_to_R(translate<RCPP_IT>(
            import_from_R(x, NA_letter),
            table,
            interpret_as_stop));
}
