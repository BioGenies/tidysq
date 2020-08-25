#ifndef TIDYSQ_COMMON_H
#define TIDYSQ_COMMON_H

#include <vector>
#include <string>
#include <Rcpp.h>

namespace tidysq::util {
    std::vector<std::string> convertStringVector(const Rcpp::StringVector &vector);
    Rcpp::StringVector convertStringVector(const std::vector<std::string> &vector);
}

#endif //TIDYSQ_COMMON_H
