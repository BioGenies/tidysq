#pragma once

#include "tidysq/tidysq-typedefs.h"

namespace tidysq::util {
    template<typename INTERNAL>
    inline LenSq random_value(LenSq alphabet_length) {
        return rand() % alphabet_length;
    }

    template<>
    inline LenSq random_value<RCPP_IT>(LenSq alphabet_length) {
        return R::runif(0, alphabet_length - 1);
    }
}