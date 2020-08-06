#ifndef TIDYSQ_OPERATIONPACK_H
#define TIDYSQ_OPERATIONPACK_H

#include "interface/Operation.h"
#include "internal/packRAWS.h"

namespace tidysq::ops {
    template<InternalType INTERNAL_IN,
            ProtoType PROTO_IN,
            InternalType INTERNAL_OUT>
    class OperationPack :
            public OperationSq<SequenceProto<INTERNAL_IN, PROTO_IN>,
                    Sequence<INTERNAL_OUT>,
                    Alphabet<INTERNAL_IN>> {
    public:
        Sequence<INTERNAL_OUT> operator()(const SequenceProto<INTERNAL_IN, PROTO_IN> &sequence,
                                          const Alphabet<INTERNAL_IN> &alphabet) const override {
            return internal::pack<INTERNAL_IN, INTERNAL_OUT>(sequence, alphabet);
        }
    };
}

#endif //TIDYSQ_OPERATIONPACK_H
