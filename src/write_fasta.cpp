#include "tidysq.h"

using namespace tidysq;

// [[Rcpp::export]]
void CPP_write_fasta(Rcpp::List &x,
                     Rcpp::StringVector &names,
                     Rcpp::StringVector &file,
                     Rcpp::IntegerVector &width,
                     Rcpp::StringVector &NA_value) {
    write_fasta(import_from_R(x, util::convert_to_scalar(NA_value)),
                util::convert_string_vector(names),
                util::convert_to_scalar(file),
                util::convert_to_scalar(width));
}