#ifndef TIDYSQ_OPERATIONUNPACK_H
#define TIDYSQ_OPERATIONUNPACK_H

#include "interface/Operation.h"
#include "internal/unpackNUMS.h"
#include "internal/unpackSTRINGS.h"
#include "internal/unpackSTRING.h"

namespace tidysq::ops {
    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT,
            ProtoType PROTO_OUT>
    class OperationUnpack;

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    class OperationUnpack<INTERNAL_IN, INTERNAL_OUT, RAWS> :
            public OperationSq<Sequence<INTERNAL_IN>,
                    SequenceProto<INTERNAL_OUT, RAWS>> {
    public:
        SequenceProto<INTERNAL_OUT, RAWS> operator()(const Sequence<INTERNAL_IN> &sequence,
                                                          const Alphabet &alphabet) const override {
            return internal::unpackNUMS<INTERNAL_IN, INTERNAL_OUT, RAWS>(sequence, alphabet);
        }
    };

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    class OperationUnpack<INTERNAL_IN, INTERNAL_OUT, INTS> :
            public OperationSq<Sequence<INTERNAL_IN>,
                    SequenceProto<INTERNAL_OUT, INTS>> {
    public:
        SequenceProto<INTERNAL_OUT, INTS> operator()(const Sequence<INTERNAL_IN> &sequence,
                                                     const Alphabet &alphabet) const override {
            return internal::unpackNUMS<INTERNAL_IN, INTERNAL_OUT, INTS>(sequence, alphabet);
        }
    };

    template<InternalType INTERNAL_IN,
            InternalType INTERNAL_OUT>
    class OperationUnpack<INTERNAL_IN, INTERNAL_OUT, STRINGS> :
            public OperationSq<Sequence<INTERNAL_IN>,
                    SequenceProto<INTERNAL_OUT, STRINGS>> {
    public:
        SequenceProto<INTERNAL_OUT, STRINGS> operator()(const Sequence<INTERNAL_IN> &sequence,
                                                        const Alphabet &alphabet) const override {
            return internal::unpackSTRINGS<INTERNAL_IN, INTERNAL_OUT>(sequence, alphabet);
        }
    };

    template<InternalType INTERNAL_IN>
    class OperationUnpack<INTERNAL_IN, STD, STRING> :
            public OperationSq<Sequence<INTERNAL_IN>,
                    SequenceProto<STD, STRING>> {
    public:
        SequenceProto<STD, STRING> operator()(const Sequence<INTERNAL_IN> &sequence,
                                              const Alphabet &alphabet) const override {
            return internal::unpackSTRING_STD<INTERNAL_IN>(sequence, alphabet);
        }
    };

    template<InternalType INTERNAL_IN>
    class OperationUnpack<INTERNAL_IN, RCPP, STRING> :
            public OperationSq<Sequence<INTERNAL_IN>,
                    SequenceProto<RCPP, STRING>> {
    public:
        SequenceProto<RCPP, STRING> operator()(const Sequence<INTERNAL_IN> &sequence,
                                               const Alphabet &alphabet) const override {
            return internal::unpackSTRING_RCPP<INTERNAL_IN>(sequence, alphabet);
        }
    };
}


#endif //TIDYSQ_OPERATIONUNPACK_H
