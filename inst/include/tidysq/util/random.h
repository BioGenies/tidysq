#pragma once

#include "tidysq/tidysq-typedefs.h"

namespace tidysq::util {
    template<typename INTERNAL>
    inline LenSq random_value(LenSq upper_bound) {
        return rand() % upper_bound;
    }

    template<>
    inline LenSq random_value<RCPP_IT>(LenSq upper_bound) {
        return R::runif(0, upper_bound - 1);
    }
}