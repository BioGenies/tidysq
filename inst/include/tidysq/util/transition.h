#ifndef TIDYSQ_TRANSITION_H
#define TIDYSQ_TRANSITION_H

#include <string>
#include <Rcpp.h>

namespace tidysq::util {
    template<int DUMMY>
    std::string getNACharacterAsString(const Rcpp::StringVector &alphabet) {
        return Rcpp::as<std::string>(Rcpp::as<Rcpp::StringVector>(alphabet.attr("na_character"))[0]);
    }
}

#endif //TIDYSQ_TRANSITION_H
