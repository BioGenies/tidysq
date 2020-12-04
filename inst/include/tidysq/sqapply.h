#pragma once

#include "tidysq/tidysq-typedefs.h"
#include "tidysq/ops/Operation.h"

namespace tidysq {
    template<typename VECTOR_IN, typename ELEMENT_IN,
             typename VECTOR_OUT, typename ELEMENT_OUT>
    VECTOR_OUT sqapply(const VECTOR_IN &vector_in,
                       ops::OperationVectorToVector<VECTOR_IN, ELEMENT_IN, VECTOR_OUT, ELEMENT_OUT> &&operation) {
        VECTOR_OUT ret = operation.initialize_vector_out(vector_in);
        for (LenSq i = 0; i < vector_in.size(); i++) {
            ret[i] = operation(vector_in[i]);
        }
        return ret;
    }

    template<typename VECTOR_IN, typename ELEMENT_IN,
            typename VECTOR_OUT, typename ELEMENT_OUT>
    VECTOR_OUT sqapply(const VECTOR_IN &vector_in,
                       ops::OperationVectorToVector<VECTOR_IN, ELEMENT_IN, VECTOR_OUT, ELEMENT_OUT> &&operation,
                       const LenSq from,
                       const LenSq to) {
        VECTOR_OUT ret = operation.initialize_vector_out(vector_in, from, to);
        for (LenSq i = 0; i < to - from; i++) {
            ret[i] = operation(vector_in[i + from]);
        }
        return ret;
    }
}
