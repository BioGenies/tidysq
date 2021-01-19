#include "tidysq/Rcpp-import.h"
#include "tidysq/Rcpp-export.h"
#include "tidysq/find_motifs.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_find_motifs(const Rcpp::List &x,
                           const std::vector<std::string> &names,
                           const std::vector<tidysq::Letter> &motifs,
                           const tidysq::Letter &NA_letter) {
    return export_to_R(find_motifs<RCPP_IT>(import_sq_from_R(x, NA_letter), names, motifs));
}
