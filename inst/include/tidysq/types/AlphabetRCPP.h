#ifndef TIDYSQ_ALPHABETRCPP_H
#define TIDYSQ_ALPHABETRCPP_H

#include <Rcpp.h>

#include "general.h"

namespace tidysq {
    template<>
    class Alphabet<RCPP> : public Rcpp::StringVector {
        using Rcpp::StringVector::StringVector;
    };
}

#endif //TIDYSQ_ALPHABETRCPP_H
