#ifndef TIDYSQ_UTILP_H
#define TIDYSQ_UTILP_H

#include "tidysq/types/Sequence.h"
#include "tidysq/types/SequenceProto.h"
#include "tidysq/types/Alphabet.h"

namespace tidysq::internal {
    template<ProtoType PROTO_OUT>
    auto matchLetter(letvalue value, const Alphabet &alphabet) -> typename SequenceProto<STD, PROTO_OUT>::ElementType;

    template<>
    inline ElemRaws matchLetter<RAWS>(const letvalue value, const Alphabet &alphabet) {
        return value;
    }

    template<InternalType INTERNAL, ProtoType PROTO>
    inline lensq calculatePackedLength(const SequenceProto<INTERNAL, PROTO> &unpacked, const Alphabet &alphabet) {
        return (alphabet.alphabetSize() * unpacked.size() + 7) / 8;
    }

    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT>
    inline Sequence<INTERNAL_OUT> reserveSpaceForPacked(const SequenceProto<INTERNAL_IN, PROTO_IN> &unpacked,
                                                        const Alphabet &alphabet) {
        return Sequence<INTERNAL_OUT>(calculatePackedLength(unpacked, alphabet), unpacked.size());
    }
}

#endif //TIDYSQ_UTILP_H
