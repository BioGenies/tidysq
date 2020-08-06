#ifndef TIDYSQ_OPERATIONUNPACK_H
#define TIDYSQ_OPERATIONUNPACK_H

#include "interface/Operation.h"
#include "internal/unpackRAWS.h"

namespace tidysq::ops {
    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT,
            ProtoType PROTO_OUT>
    class OperationUnpack :
            public OperationSq<Sequence<INTERNAL_IN>,
                    SequenceProto<INTERNAL_OUT, PROTO_OUT>,
                    Alphabet<INTERNAL_IN>> {
    public:
        SequenceProto<INTERNAL_OUT, PROTO_OUT> operator()(const Sequence<INTERNAL_IN> &sequence,
                                                          const Alphabet<INTERNAL_IN> &alphabet) const override {
            return internal::unpack<INTERNAL_IN, INTERNAL_OUT>(sequence, alphabet);
        }
    };
}


#endif //TIDYSQ_OPERATIONUNPACK_H
