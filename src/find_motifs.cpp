#include <Rcpp.h>

#include "tidysq/exports.h"
#include "tidysq/ops/find_motifs.h"

using namespace tidysq;

//[[Rcpp::export]]
Rcpp::List CPP_find_motifs(const Rcpp::List &x,
                           const Rcpp::StringVector &names,
                           const Rcpp::StringVector &motifs) {
    return find_motifs<RCPP>(importFromR(x, "!"),
                             Rcpp::as<std::vector<std::string>>(names),
                             Rcpp::as<std::vector<std::string>>(motifs)).exportToR();
}

//[[Rcpp::export]]
Rcpp::LogicalVector CPP_has(const Rcpp::List &x,
                            const Rcpp::StringVector &motifs) {
    return has<RCPP>(importFromR(x, "!"), Rcpp::as<std::vector<std::string>>(motifs));
}
