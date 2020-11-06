#pragma once

#include "tidysq/types/Sequence.h"
#include "tidysq/types/ProtoSequence.h"
#include "tidysq/types/Alphabet.h"

namespace tidysq::internal {
    template<ProtoType PROTO_OUT>
    auto match_letter(LetterValue value, const Alphabet &alphabet)
            -> typename ProtoTypeMapper<PROTO_OUT>::ProtoSequenceElementType;

    template<>
    inline auto match_letter<RAWS>(const LetterValue value, const Alphabet &alphabet)
            -> typename ProtoTypeMapper<RAWS>::ProtoSequenceElementType {
        return value;
    }

    template<>
    inline auto match_letter<INTS>(const LetterValue value, const Alphabet &alphabet)
            -> typename ProtoTypeMapper<INTS>::ProtoSequenceElementType {
        return value;
    }

    template<>
    inline auto match_letter<STRINGS>(const LetterValue value, const Alphabet &alphabet)
            -> typename ProtoTypeMapper<STRINGS>::ProtoSequenceElementType {
        return alphabet[value];
    }

    template<>
    inline auto match_letter<STRING>(const LetterValue value, const Alphabet &alphabet)
            -> typename ProtoTypeMapper<STRING>::ProtoSequenceElementType {
        return alphabet.get_simple_letter(value);
    }

    inline Letter match_letter_multichar(const LetterValue value, const Alphabet &alphabet) {
        return alphabet[value];
    }

    inline LenSq calculate_packed_internal_length(const LenSq unpackedLength, const AlphSize &alph_size) {
        return (alph_size * unpackedLength + 7) / 8;
    }

    template<InternalType INTERNAL_OUT>
    inline Sequence<INTERNAL_OUT> reserve_space_for_packed(const LenSq &unpacked_length,
                                                           const AlphSize &alph_size) {
        return Sequence<INTERNAL_OUT>(calculate_packed_internal_length(unpacked_length, alph_size), unpacked_length);
    }

    template<InternalType INTERNAL_IN, InternalType INTERNAL_OUT, ProtoType PROTO_OUT>
    inline ProtoSequence<INTERNAL_OUT, PROTO_OUT> reserve_space_for_unpacked(const Sequence<INTERNAL_IN> &packed) {
        return ProtoSequence<INTERNAL_OUT, PROTO_OUT>(packed.originalLength());
    }
}
