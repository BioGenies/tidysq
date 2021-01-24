#pragma once

#include "tidysq/tidysq-typedefs.h"

namespace tidysq::util {
    template<typename INTERNAL>
    inline LenSq find_common_size(const std::vector<Sq<INTERNAL>> &list_of_sq) {
        LenSq common_size = 1;
        for (const Sq<INTERNAL>& sq : list_of_sq) {
            if (sq.size() != 1) {
                if (common_size == 1) {
                    common_size = sq.size();
                } else if (common_size != sq.size()) {
                    throw std::invalid_argument("all sq objects must have either the same size or size 1");
                }
            }
        }
        return common_size;
    }
}
