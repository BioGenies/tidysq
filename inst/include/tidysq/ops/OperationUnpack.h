#ifndef TIDYSQ_OPERATIONUNPACK_H
#define TIDYSQ_OPERATIONUNPACK_H

#include "interface/Operation.h"
#include "internal/unpackNUMS.h"
#include "internal/unpackSTRINGS.h"

namespace tidysq::ops {
    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT,
            ProtoType PROTO_OUT>
    class OperationUnpack;

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    class OperationUnpack<INTERNAL_IN, INTERNAL_OUT, RAWS> :
            public OperationSq<Sequence<INTERNAL_IN>,
                    SequenceProto<INTERNAL_OUT, RAWS>,
                    Alphabet<INTERNAL_IN>> {
    public:
        SequenceProto<INTERNAL_OUT, RAWS> operator()(const Sequence<INTERNAL_IN> &sequence,
                                                          const Alphabet<INTERNAL_IN> &alphabet) const override {
            return internal::unpackNUMS<INTERNAL_IN, INTERNAL_OUT, RAWS>(sequence, alphabet);
        }
    };

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    class OperationUnpack<INTERNAL_IN, INTERNAL_OUT, INTS> :
            public OperationSq<Sequence<INTERNAL_IN>,
                    SequenceProto<INTERNAL_OUT, INTS>,
                    Alphabet<INTERNAL_IN>> {
    public:
        SequenceProto<INTERNAL_OUT, INTS> operator()(const Sequence<INTERNAL_IN> &sequence,
                                                     const Alphabet<INTERNAL_IN> &alphabet) const override {
            return internal::unpackNUMS<INTERNAL_IN, INTERNAL_OUT, INTS>(sequence, alphabet);
        }
    };

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    class OperationUnpack<INTERNAL_IN, INTERNAL_OUT, STRINGS> :
            public OperationSq<Sequence<INTERNAL_IN>,
                    SequenceProto<INTERNAL_OUT, STRINGS>,
                    Alphabet<INTERNAL_IN>> {
    public:
        SequenceProto<INTERNAL_OUT, STRINGS> operator()(const Sequence<INTERNAL_IN> &sequence,
                                                        const Alphabet<INTERNAL_IN> &alphabet) const override {
            return internal::unpackSTRINGS<INTERNAL_IN, INTERNAL_OUT>(sequence, alphabet);
        }
    };
}


#endif //TIDYSQ_OPERATIONUNPACK_H
