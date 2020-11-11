#ifndef TIDYSQ_OPERATIONUNPACK_H
#define TIDYSQ_OPERATIONUNPACK_H

#include "tidysq/ops/Operation.h"
#include "tidysq/internal/unpack_common.h"
#include "tidysq/internal/unpack_multichar.h"
#include "tidysq/util/calculate_length.h"

namespace tidysq::ops {
    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT, ProtoType PROTO_OUT>
    class OperationUnpack : public tidysq::ops::OperationSq<Sequence<INTERNAL_IN>,
                                               ProtoSequence<INTERNAL_OUT, PROTO_OUT>> {
    public:
        ProtoSequence<INTERNAL_OUT, PROTO_OUT> operator()(const Sequence<INTERNAL_IN> &packed,
                                                          const Alphabet &alphabet) const override {
            ProtoSequence<INTERNAL_OUT, PROTO_OUT> unpacked = util::reserve_space_for_unpacked<INTERNAL_IN, INTERNAL_OUT, PROTO_OUT>(
                    packed);
            internal::unpack_common<INTERNAL_IN, INTERNAL_OUT, PROTO_OUT>(packed, unpacked, alphabet);
            return unpacked;
        }
    };

    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT>
    class OperationUnpack<INTERNAL_IN, INTERNAL_OUT, STRING> : public tidysq::ops::OperationSq<Sequence<INTERNAL_IN>,
            ProtoSequence<INTERNAL_OUT, STRING>> {
    public:
        ProtoSequence<INTERNAL_OUT, STRING> operator()(const Sequence<INTERNAL_IN> &packed,
                                                          const Alphabet &alphabet) const override {
            ProtoSequence<INTERNAL_OUT, STRING> unpacked;
            if (alphabet.is_simple()) {
                unpacked = util::reserve_space_for_unpacked<INTERNAL_IN, INTERNAL_OUT, STRING>(packed);
                internal::unpack_common<INTERNAL_IN, INTERNAL_OUT, STRING>(packed, unpacked, alphabet);
            } else {
                unpacked = {};
                internal::unpack_multichar_string<INTERNAL_IN, INTERNAL_OUT>(packed, unpacked, alphabet);
            }
            return unpacked;
        }
    };
}


#endif //TIDYSQ_OPERATIONUNPACK_H
