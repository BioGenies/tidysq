#ifndef TIDYSQ_UTIL_H
#define TIDYSQ_UTIL_H

#include "tidysq/types/Sequence.h"
#include "tidysq/types/ProtoSequence.h"
#include "tidysq/types/Alphabet.h"

namespace tidysq::internal {
    template<ProtoType PROTO_OUT>
    auto matchLetter(LetterValue value, const Alphabet &alphabet)
            -> typename ProtoTypeMapper<PROTO_OUT>::ProtoSequenceElementType;

    template<>
    inline auto matchLetter<RAWS>(const LetterValue value, const Alphabet &alphabet)
            -> typename ProtoTypeMapper<RAWS>::ProtoSequenceElementType {
        return value;
    }

    template<>
    inline auto matchLetter<INTS>(const LetterValue value, const Alphabet &alphabet)
            -> typename ProtoTypeMapper<INTS>::ProtoSequenceElementType {
        return value;
    }

    template<>
    inline auto matchLetter<STRINGS>(const LetterValue value, const Alphabet &alphabet)
            -> typename ProtoTypeMapper<STRINGS>::ProtoSequenceElementType {
        return alphabet[value];
    }

    template<>
    inline auto matchLetter<STRING>(const LetterValue value, const Alphabet &alphabet)
            -> typename ProtoTypeMapper<STRING>::ProtoSequenceElementType {
        return alphabet.get_simple_letter(value);
    }

    inline Letter matchLetterMultichar(const LetterValue value, const Alphabet &alphabet) {
        return alphabet[value];
    }

    inline LenSq calculate_packed_internal_length(const LenSq unpackedLength, const Alphabet &alphabet) {
        return (alphabet.alphabet_size() * unpackedLength + 7);
    }

    template<InternalType INTERNAL, ProtoType PROTO>
    inline LenSq calculate_packed_internal_length(const ProtoSequence<INTERNAL, PROTO> &unpacked, const Alphabet &alphabet) {
        return calculate_packed_internal_length(unpacked.length(), alphabet);
    }

    template<InternalType INTERNAL_IN, ProtoType PROTO_IN, InternalType INTERNAL_OUT>
    inline Sequence<INTERNAL_OUT> reserveSpaceForPacked(const ProtoSequence<INTERNAL_IN, PROTO_IN> &unpacked,
                                                        const Alphabet &alphabet) {
        return Sequence<INTERNAL_OUT>(calculate_packed_internal_length(unpacked, alphabet), unpacked.length());
    }

    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT, ProtoType PROTO_OUT>
    inline ProtoSequence<INTERNAL_OUT, PROTO_OUT> reserveSpaceForUnpacked(const Sequence<INTERNAL_IN> &packed) {
        return ProtoSequence<INTERNAL_OUT, PROTO_OUT>(packed.originalLength());
    }
}

#endif //TIDYSQ_UTIL_H
