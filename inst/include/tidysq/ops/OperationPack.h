#ifndef TIDYSQ_OPERATIONPACK_H
#define TIDYSQ_OPERATIONPACK_H

#include "tidysq/ops/interface/Operation.h"
#include "tidysq/ops/internal/pack.h"
#include "tidysq/ops/internal/util.h"

namespace tidysq::ops {
    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT>
    class OperationPack : public OperationSq<ProtoSequence<INTERNAL_IN, PROTO_IN>,
                                             Sequence<INTERNAL_OUT>> {
    public:
        Sequence<INTERNAL_OUT> operator() (const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
                                           const Alphabet &alphabet) const override {
            Sequence<INTERNAL_OUT> packed = internal::reserveSpaceForPacked<INTERNAL_IN, PROTO_IN, INTERNAL_OUT>(unpacked, alphabet);
            if (alphabet.is_simple()) {
                internal::pack<INTERNAL_IN, PROTO_IN, INTERNAL_OUT, true>(unpacked, packed, alphabet);
            } else {
                internal::pack<INTERNAL_IN, PROTO_IN, INTERNAL_OUT, false>(unpacked, packed, alphabet);
            }
            return packed;
        }
    };
}

#endif //TIDYSQ_OPERATIONPACK_H
