#include <Rcpp.h>

#include "tidysq/exports.h"

using namespace tidysq;

//' @export
//[[Rcpp::export]]
Rcpp::DataFrame CPP_read_fasta(Rcpp::StringVector file_name,
                          Rcpp::StringVector alphabet) {
  auto val = readFasta<RCPP>(util::getScalarStringValue(file_name), Alphabet(alphabet));
  auto ret = Rcpp::DataFrame::create(Rcpp::Named("sq") = std::get<0>(val).exportToR(),
                                     Rcpp::Named("name") = util::convertStringVector(std::get<1>(val)));
  ret.attr("class") = Rcpp::StringVector{"tbl_df", "tbl", "data.frame"};
  return ret;
}