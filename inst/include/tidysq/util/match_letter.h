#pragma once

#include "tidysq/Sequence.h"
#include "tidysq/ProtoSequence.h"
#include "tidysq/Alphabet.h"

namespace tidysq::util {
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
}

