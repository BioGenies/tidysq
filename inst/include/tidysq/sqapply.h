#pragma once

#include "tidysq/tidysq-typedefs.h"
#include "tidysq/ops/Operation.h"

namespace tidysq {
    template <typename TYPE_IN, typename TYPE_OUT>
    TYPE_OUT sqapply(const TYPE_IN &sq,
                     const ops::OperationSq<typename TYPE_IN::ElementType, typename TYPE_OUT::ElementType> &op) {
        TYPE_OUT ret(sq.length(), sq.alphabet());
        for (LenSq i = 0; i < sq.length(); i++) {
            ret[i] = op(sq[i], sq.alphabet());
        }
        return ret;
    }

    template <typename TYPE_IN, typename TYPE_OUT>
    TYPE_OUT sqapply(const TYPE_IN &sq,
                     const ops::OperationSq<typename TYPE_IN::ElementType, typename TYPE_OUT::ElementType> &op,
                     const LenSq from,
                     const LenSq to) {
        TYPE_OUT ret(to - from, sq.alphabet());
        for (LenSq i = 0; i < to - from; i++) {
            ret[i] = op(sq[from + i], sq.alphabet());
        }
        return ret;
    }
}
