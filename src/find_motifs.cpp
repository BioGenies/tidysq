#include "tidysq/Rcpp-import.h"
#include "tidysq/Rcpp-export.h"
#include "tidysq/find_motifs.h"

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

