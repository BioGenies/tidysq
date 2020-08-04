#ifndef TIDYSQ_IMPL_OPERATION_H
#define TIDYSQ_IMPL_OPERATION_H

#include "interface_Operation.h"
#include "internal_pack.h"

namespace tidysq {
    namespace ops {
        template <typename SEQUENCE_IN, typename SEQUENCE_OUT, typename ALPHABET>
        class OperationPack : public OperationSq<SEQUENCE_IN, SEQUENCE_OUT, ALPHABET> {
        public:
            SEQUENCE_OUT operator() (const SEQUENCE_IN &sequence, const ALPHABET &alphabet) const override {
                return internal::pack<SEQUENCE_IN, SEQUENCE_OUT, ALPHABET>(sequence, alphabet);
            }
        };
    }
}

#endif //TIDYSQ_IMPL_OPERATION_H
