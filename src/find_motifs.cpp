#include <Rcpp.h>

#include "tidysq.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_find_motifs(const Rcpp::List &x,
                           const Rcpp::StringVector &names,
                           const Rcpp::StringVector &motifs,
                           const Rcpp::StringVector &NA_letter) {
    return export_to_R(find_motifs<RCPP_IT>(import_from_R(x, NA_letter),
                                         util::convert_string_vector(names),
                                         util::convert_string_vector(motifs)));
}

//[[Rcpp::export]]
Rcpp::LogicalVector CPP_has(const Rcpp::List &x,
                            const Rcpp::StringVector &motifs,
                            const Rcpp::StringVector &NA_letter) {
   return has<RCPP_IT>(import_from_R(x, NA_letter), util::convert_string_vector(motifs));
}
