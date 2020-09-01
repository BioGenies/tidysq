#ifndef TIDYSQ_SQAPPLY_H
#define TIDYSQ_SQAPPLY_H

#include "types/general.h"
#include "ops/interface/Operation.h"

namespace tidysq {
    template <typename TYPE_IN, typename TYPE_OUT>
    TYPE_OUT sqapply(const TYPE_IN &sq, const OperationSq<typename TYPE_IN::SequenceType, typename TYPE_OUT::SequenceType> &op) {
        TYPE_OUT ret(sq.length(), sq.alphabet(), sq.type());
        for (lensq i = 0; i < sq.length(); i++) {
            ret[i] = op(sq[i], sq.alphabet());
        }
        return ret;
    }
}

#endif //TIDYSQ_SQAPPLY_H
