#ifndef TIDYSQ_OPERATIONPACK_H
#define TIDYSQ_OPERATIONPACK_H

#include "interface/Operation.h"
#include "internal/packNUMS.h"
#include "internal/packSTRINGS.h"
#include "internal/packSTRING.h"

namespace tidysq::ops {
    template<InternalType INTERNAL_IN,
            ProtoType PROTO_IN,
            InternalType INTERNAL_OUT>
    class OperationPack;

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    class OperationPack<INTERNAL_IN, RAWS, INTERNAL_OUT> :
            public OperationSq<SequenceProto<INTERNAL_IN, RAWS>,
                    Sequence<INTERNAL_OUT>> {
    public:
        Sequence<INTERNAL_OUT> operator()(const SequenceProto<INTERNAL_IN, RAWS> &sequence,
                                          const Alphabet &alphabet) const override {
            return internal::packNUMS<INTERNAL_IN, RAWS, INTERNAL_OUT>(sequence, alphabet);
        }
    };

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    class OperationPack<INTERNAL_IN, INTS, INTERNAL_OUT> :
            public OperationSq<SequenceProto<INTERNAL_IN, INTS>,
                    Sequence<INTERNAL_OUT>> {
    public:
        Sequence<INTERNAL_OUT> operator()(const SequenceProto<INTERNAL_IN, INTS> &sequence,
                                          const Alphabet &alphabet) const override {
            return internal::packNUMS<INTERNAL_IN, INTS, INTERNAL_OUT>(sequence, alphabet);
        }
    };

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    class OperationPack<INTERNAL_IN, STRINGS, INTERNAL_OUT> :
            public OperationSq<SequenceProto<INTERNAL_IN, STRINGS>,
                    Sequence<INTERNAL_OUT>> {
    public:
        Sequence<INTERNAL_OUT> operator()(const SequenceProto<INTERNAL_IN, STRINGS> &sequence,
                                          const Alphabet &alphabet) const override {
            return internal::packSTRINGS<INTERNAL_IN, INTERNAL_OUT>(sequence, alphabet);
        }
    };

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    class OperationPack<INTERNAL_IN, STRING, INTERNAL_OUT> :
            public OperationSq<SequenceProto<INTERNAL_IN, STRING>,
                    Sequence<INTERNAL_OUT>> {
    public:
        Sequence<INTERNAL_OUT> operator()(const SequenceProto<INTERNAL_IN, STRING> &sequence,
                                          const Alphabet &alphabet) const override {
            return internal::packSTRING<INTERNAL_IN, INTERNAL_OUT>(sequence, alphabet);
        }
    };
}

#endif //TIDYSQ_OPERATIONPACK_H
