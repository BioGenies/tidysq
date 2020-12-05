#include "tidysq.h"

using namespace tidysq;

// [[Rcpp::export]]
void CPP_write_fasta(const Rcpp::List &x,
                     const std::vector<std::string> &names,
                     const std::string &file,
                     const int &width,
                     const tidysq::Letter &NA_value) {
    write_fasta(import_from_R(x, NA_value), names, file, width);
}