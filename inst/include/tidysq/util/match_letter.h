#pragma once

#include "tidysq/Sequence.h"
#include "tidysq/ProtoSequence.h"
#include "tidysq/Alphabet.h"

namespace tidysq::util {
    template<typename PROTO_OUT>
    auto match_letter(LetterValue value, [[maybe_unused]] const Alphabet &alphabet)
            -> typename ProtoTypeMapper<PROTO_OUT>::ProtoSequenceElementType;

    template<>
    inline auto match_letter<RAWS_PT>(const LetterValue value, [[maybe_unused]] const Alphabet &alphabet)
            -> typename ProtoTypeMapper<RAWS_PT>::ProtoSequenceElementType {
        return value;
    }

    template<>
    inline auto match_letter<INTS_PT>(const LetterValue value, [[maybe_unused]] const Alphabet &alphabet)
            -> typename ProtoTypeMapper<INTS_PT>::ProtoSequenceElementType {
        return value;
    }

    template<>
    inline auto match_letter<STRINGS_PT>(const LetterValue value, const Alphabet &alphabet)
            -> typename ProtoTypeMapper<STRINGS_PT>::ProtoSequenceElementType {
        return alphabet[value];
    }

    template<>
    inline auto match_letter<STRING_PT>(const LetterValue value, const Alphabet &alphabet)
            -> typename ProtoTypeMapper<STRING_PT>::ProtoSequenceElementType {
        return alphabet.get_simple_letter(value);
    }

    inline Letter match_letter_multichar(const LetterValue value, const Alphabet &alphabet) {
        return alphabet[value];
    }
}

