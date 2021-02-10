#pragma once

#include "tidysq/tidysq-typedefs.h"

namespace tidysq::constants {
    template<bool SIMPLE, LetterValue CODON>
    const LetterValue COMPLEMENT = CODON;

    template<bool SIMPLE>
    const LetterValue COMPLEMENT<SIMPLE, 0u> = 3u;

    template<bool SIMPLE>
    const LetterValue COMPLEMENT<SIMPLE, 1u> = 2u;

    template<bool SIMPLE>
    const LetterValue COMPLEMENT<SIMPLE, 2u> = 1u;

    template<bool SIMPLE>
    const LetterValue COMPLEMENT<SIMPLE, 3u> = 0u;

    template<>
    const LetterValue COMPLEMENT<false, 6u> = 7u;

    template<>
    const LetterValue COMPLEMENT<false, 7u> = 6u;

    template<>
    const LetterValue COMPLEMENT<false, 8u> = 9u;

    template<>
    const LetterValue COMPLEMENT<false, 9u> = 8u;

    template<>
    const LetterValue COMPLEMENT<false, 10u> = 13u;

    template<>
    const LetterValue COMPLEMENT<false, 11u> = 12u;

    template<>
    const LetterValue COMPLEMENT<false, 12u> = 11u;

    template<>
    const LetterValue COMPLEMENT<false, 13u> = 10u;
}
