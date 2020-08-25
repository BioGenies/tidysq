#include "tidysq/util/common.h"

namespace tidysq::util {
std::vector<std::string> convertStringVector(const Rcpp::StringVector &vector) {
  std::vector<std::string> ret(vector.size());
  auto iterIn = vector.begin();
  auto iterOut = ret.begin();
  while (iterIn != vector.end()) {
    *iterOut = *iterIn;
    iterIn++;
    iterOut++;
  }
  return ret;
}

Rcpp::StringVector convertStringVector(const std::vector<std::string> &vector) {
  Rcpp::StringVector ret(vector.size());
  auto iterIn = vector.begin();
  auto iterOut = ret.begin();
  while (iterIn != vector.end()) {
    *iterOut = *iterIn;
    iterIn++;
    iterOut++;
  }
  return ret;
}
}