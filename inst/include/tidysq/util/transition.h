#ifndef TIDYSQ_TRANSITION_H
#define TIDYSQ_TRANSITION_H

#include <string>
#include <Rcpp.h>

namespace tidysq::util {
    std::string getNACharacterAsString(const Rcpp::StringVector &alphabet);
}

#endif //TIDYSQ_TRANSITION_H
