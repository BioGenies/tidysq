#ifndef TIDYSQ_UTIL_H
#define TIDYSQ_UTIL_H

#include "tidysq/types/Sequence.h"
#include "tidysq/types/ProtoSequence.h"
#include "tidysq/types/Alphabet.h"

namespace tidysq::internal {
    template<ProtoType PROTO_OUT>
    auto matchLetter(LetValue value, const Alphabet &alphabet) -> typename TypeMapper<STD, PROTO_OUT>::ProtoSequenceElementType;

    template<>
    inline ElemRaws matchLetter<RAWS>(const LetValue value, const Alphabet &alphabet) {
        return value;
    }

    template<>
    inline ElemInts matchLetter<INTS>(const LetValue value, const Alphabet &alphabet) {
        return value;
    }

    template<>
    inline ElemStrings matchLetter<STRINGS>(const LetValue value, const Alphabet &alphabet) {
        return alphabet[value];
    }

    template<InternalType INTERNAL, ProtoType PROTO>
    inline LenSq calculatePackedLength(const ProtoSequence<INTERNAL, PROTO> &unpacked, const Alphabet &alphabet) {
        return (alphabet.alphabetSize() * unpacked.length() + 7) / 8;
    }

    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT>
    inline Sequence<INTERNAL_OUT> reserveSpaceForPacked(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
                                                        const Alphabet &alphabet) {
        return Sequence<INTERNAL_OUT>(calculatePackedLength(unpacked, alphabet), unpacked.length());
    }

    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT, ProtoType PROTO_OUT>
    inline ProtoSequence<INTERNAL_OUT, PROTO_OUT> reserveSpaceForUnpacked(const Sequence<INTERNAL_IN> &packed) {
        return ProtoSequence<INTERNAL_OUT, PROTO_OUT>(packed.originalLength());
    }

    template<InternalType INTERNAL, ProtoType PROTO>
    inline LetValue matchValue(const typename TypeMapper<INTERNAL, PROTO>::ProtoSequenceElementType &element, const Alphabet &alphabet) {
        for (unsigned int i = 0; i < alphabet.length(); i++) {
            if (element == alphabet[i]) return i;
        }
        return alphabet.NAValue();
    }
}

#endif //TIDYSQ_UTIL_H
