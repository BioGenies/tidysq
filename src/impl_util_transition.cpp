#include "tidysq/util/transition.h"

namespace tidysq::util {
  std::string getNACharacterAsString(const Rcpp::StringVector &alphabet) {
    return Rcpp::as<std::string>(Rcpp::as<Rcpp::StringVector>(alphabet.attr("na_character"))[0]);
  }
}