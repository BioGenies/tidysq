#include <Rcpp.h>

#include "tidysq.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_find_motifs(const Rcpp::List &x,
                           const std::vector<std::string> &names,
                           const std::vector<tidysq::Letter> &motifs,
                           const tidysq::Letter &NA_letter) {
    return export_to_R(find_motifs<RCPP_IT>(import_from_R(x, NA_letter), names, motifs));
}

//[[Rcpp::export]]
Rcpp::LogicalVector CPP_has(const Rcpp::List &x,
                            const std::vector<tidysq::Letter> &motifs,
                            const tidysq::Letter &NA_letter) {
   return has<RCPP_IT>(import_from_R(x, NA_letter), motifs);
}
