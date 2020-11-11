#ifndef TIDYSQ_OPERATIONPACK_H
#define TIDYSQ_OPERATIONPACK_H

#include "tidysq/ops/Operation.h"
#include "tidysq/internal/pack.h"
#include "tidysq/util/calculate_length.h"

namespace tidysq::ops {
    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT>
    class OperationPack : public tidysq::ops::OperationSq<ProtoSequence<INTERNAL_IN, PROTO_IN>,
                                             Sequence<INTERNAL_OUT>> {
    public:
        Sequence<INTERNAL_OUT> operator() (const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
                                           const Alphabet &alphabet) const override {
            Sequence<INTERNAL_OUT> packed = util::reserve_space_for_packed<INTERNAL_OUT>(
                    unpacked.length(), alphabet.alphabet_size());
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
