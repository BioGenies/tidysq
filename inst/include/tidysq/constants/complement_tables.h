#pragma once

#include "tidysq/tidysq-typedefs.h"

namespace tidysq::constants {
    template<bool SIMPLE, LetterValue CODON>
    inline const LetterValue COMPLEMENT = CODON;

    template<bool SIMPLE>
    inline const LetterValue COMPLEMENT<SIMPLE, 0u> = 3u;

    template<bool SIMPLE>
    inline const LetterValue COMPLEMENT<SIMPLE, 1u> = 2u;

    template<bool SIMPLE>
    inline const LetterValue COMPLEMENT<SIMPLE, 2u> = 1u;

    template<bool SIMPLE>
    inline const LetterValue COMPLEMENT<SIMPLE, 3u> = 0u;

    template<>
    inline const LetterValue COMPLEMENT<false, 6u> = 7u;

    template<>
    inline const LetterValue COMPLEMENT<false, 7u> = 6u;

    template<>
    inline const LetterValue COMPLEMENT<false, 8u> = 9u;

    template<>
    inline const LetterValue COMPLEMENT<false, 9u> = 8u;

    template<>
    inline const LetterValue COMPLEMENT<false, 10u> = 13u;

    template<>
    inline const LetterValue COMPLEMENT<false, 11u> = 12u;

    template<>
    inline const LetterValue COMPLEMENT<false, 12u> = 11u;

    template<>
    inline const LetterValue COMPLEMENT<false, 13u> = 10u;
}
