#ifndef TIDYSQ_OPERATIONUNPACK_H
#define TIDYSQ_OPERATIONUNPACK_H

#include "tidysq/ops/interface/Operation.h"
#include "tidysq/ops/internal/unpackStandard.h"
#include "tidysq/ops/internal/util.h"

namespace tidysq::ops {
    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT, ProtoType PROTO_OUT>
    class OperationUnpack : public OperationSq<Sequence<INTERNAL_IN>,
                                               ProtoSequence<INTERNAL_OUT, PROTO_OUT>> {
    public:
        ProtoSequence<INTERNAL_OUT, PROTO_OUT> operator()(const Sequence<INTERNAL_IN> &packed,
                                                          const Alphabet &alphabet) const override {
            ProtoSequence<INTERNAL_OUT, PROTO_OUT> unpacked = internal::reserveSpaceForUnpacked<INTERNAL_IN, INTERNAL_OUT, PROTO_OUT>(packed);
            internal::unpackStandard<INTERNAL_IN, INTERNAL_OUT, PROTO_OUT>(packed, unpacked, alphabet);
            return unpacked;
        }
    };
}


#endif //TIDYSQ_OPERATIONUNPACK_H
