#pragma once

#include "tidysq/tidysq-typedefs.h"

namespace tidysq::constants {
    template<int TABLE, LetterValue C1, LetterValue C2, LetterValue C3>
    inline const LetterValue CODON = 31u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 0u, 0u, 0u> = 8u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 0u, 0u, 1u> = 11u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 0u, 0u, 2u> = 8u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 0u, 0u, 3u> = 11u;

    template<int TABLE, LetterValue C3>
    inline const LetterValue CODON<TABLE, 0u, 1u, C3> = 16u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 0u, 2u, 0u> = 14u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 0u, 2u, 1u> = 15u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 0u, 2u, 2u> = 14u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 0u, 2u, 3u> = 15u;

    template<int TABLE, LetterValue C3>
    inline const LetterValue CODON<TABLE, 0u, 3u, C3> = 7u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 0u, 3u, 2u> = 10u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 1u, 0u, 0u> = 13u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 1u, 0u, 1u> = 6u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 1u, 0u, 2u> = 13u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 1u, 0u, 3u> = 6u;

    template<int TABLE, LetterValue C3>
    inline const LetterValue CODON<TABLE, 1u, 1u, C3> = 12u;

    template<int TABLE, LetterValue C3>
    inline const LetterValue CODON<TABLE, 1u, 2u, C3> = 14u;

    template<int TABLE, LetterValue C3>
    inline const LetterValue CODON<TABLE, 1u, 3u, C3> = 9u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 2u, 0u, 0u> = 3u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 2u, 0u, 1u> = 2u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 2u, 0u, 2u> = 3u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 2u, 0u, 3u> = 2u;

    template<int TABLE, LetterValue C3>
    inline const LetterValue CODON<TABLE, 2u, 1u, C3> = 0u;

    template<int TABLE, LetterValue C3>
    inline const LetterValue CODON<TABLE, 2u, 2u, C3> = 5u;

    template<int TABLE, LetterValue C3>
    inline const LetterValue CODON<TABLE, 2u, 3u, C3> = 17u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 3u, 0u, 0u> = 21u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 3u, 0u, 1u> = 19u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 3u, 0u, 2u> = 21u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 3u, 0u, 3u> = 19u;

    template<int TABLE, LetterValue C3>
    inline const LetterValue CODON<TABLE, 3u, 1u, C3> = 15u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 3u, 2u, 0u> = 21u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 3u, 2u, 1u> = 1u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 3u, 2u, 2u> = 18u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 3u, 2u, 3u> = 1u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 3u, 3u, 0u> = 4u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 3u, 3u, 1u> = 9u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 3u, 3u, 2u> = 4u;

    template<int TABLE>
    inline const LetterValue CODON<TABLE, 3u, 3u, 3u> = 9u;


    template<>
    inline const LetterValue CODON<2, 0u, 2u, 0u> = 21u;

    template<>
    inline const LetterValue CODON<2, 0u, 2u, 2u> = 21u;

    template<>
    inline const LetterValue CODON<2, 0u, 3u, 0u> = 10u;

    template<>
    inline const LetterValue CODON<2, 3u, 2u, 0u> = 18u;


    template<>
    inline const LetterValue CODON<3, 0u, 3u, 0u> = 10u;

    template<>
    inline const LetterValue CODON<3, 1u, 2u, 0u> = 31u;

    template<>
    inline const LetterValue CODON<3, 1u, 2u, 1u> = 31u;

    template<LetterValue C3>
    inline const LetterValue CODON<3, 1u, 3u, C3> = 16u;

    template<>
    inline const LetterValue CODON<3, 3u, 2u, 0u> = 18u;


    template<>
    inline const LetterValue CODON<4, 3u, 2u, 0u> = 18u;


    template<>
    inline const LetterValue CODON<5, 0u, 2u, 0u> = 15u;

    template<>
    inline const LetterValue CODON<5, 0u, 2u, 2u> = 15u;

    template<>
    inline const LetterValue CODON<5, 0u, 3u, 0u> = 10u;

    template<>
    inline const LetterValue CODON<5, 3u, 2u, 0u> = 18u;


    template<>
    inline const LetterValue CODON<6, 3u, 0u, 0u> = 13u;

    template<>
    inline const LetterValue CODON<6, 3u, 0u, 2u> = 13u;


    template<>
    inline const LetterValue CODON<9, 0u, 0u, 0u> = 11u;

    template<>
    inline const LetterValue CODON<9, 0u, 2u, 0u> = 15u;

    template<>
    inline const LetterValue CODON<9, 0u, 2u, 2u> = 15u;

    template<>
    inline const LetterValue CODON<9, 3u, 2u, 0u> = 18u;


    template<>
    inline const LetterValue CODON<10, 3u, 2u, 0u> = 1u;


    template<>
    inline const LetterValue CODON<12, 1u, 3u, 2u> = 15u;


    template<>
    inline const LetterValue CODON<13, 0u, 2u, 0u> = 5u;

    template<>
    inline const LetterValue CODON<13, 0u, 2u, 2u> = 5u;

    template<>
    inline const LetterValue CODON<13, 0u, 3u, 0u> = 10u;

    template<>
    inline const LetterValue CODON<13, 3u, 2u, 0u> = 18u;


    template<>
    inline const LetterValue CODON<14, 0u, 0u, 0u> = 11u;

    template<>
    inline const LetterValue CODON<14, 0u, 2u, 0u> = 15u;

    template<>
    inline const LetterValue CODON<14, 0u, 2u, 2u> = 15u;

    template<>
    inline const LetterValue CODON<14, 3u, 0u, 0u> = 19u;

    template<>
    inline const LetterValue CODON<14, 3u, 2u, 0u> = 18u;


    template<>
    inline const LetterValue CODON<15, 3u, 0u, 2u> = 13u;


    template<>
    inline const LetterValue CODON<16, 3u, 0u, 2u> = 9u;


    template<>
    inline const LetterValue CODON<21, 0u, 0u, 0u> = 11u;

    template<>
    inline const LetterValue CODON<21, 0u, 2u, 0u> = 15u;

    template<>
    inline const LetterValue CODON<21, 0u, 2u, 2u> = 15u;

    template<>
    inline const LetterValue CODON<21, 0u, 3u, 0u> = 10u;

    template<>
    inline const LetterValue CODON<21, 3u, 2u, 0u> = 18u;


    template<>
    inline const LetterValue CODON<22, 3u, 0u, 2u> = 9u;

    template<>
    inline const LetterValue CODON<22, 3u, 1u, 0u> = 21u;


    template<>
    inline const LetterValue CODON<23, 3u, 3u, 0u> = 21u;


    template<>
    inline const LetterValue CODON<24, 0u, 2u, 0u> = 15u;

    template<>
    inline const LetterValue CODON<24, 0u, 2u, 2u> = 8u;

    template<>
    inline const LetterValue CODON<24, 3u, 2u, 0u> = 18u;


    template<>
    inline const LetterValue CODON<25, 3u, 2u, 0u> = 5u;


    template<>
    inline const LetterValue CODON<26, 1u, 3u, 2u> = 0u;


    template<>
    inline const LetterValue CODON<29, 3u, 0u, 0u> = 19u;

    template<>
    inline const LetterValue CODON<29, 3u, 0u, 2u> = 19u;


    template<>
    inline const LetterValue CODON<30, 3u, 0u, 0u> = 3u;

    template<>
    inline const LetterValue CODON<30, 3u, 0u, 2u> = 3u;


    template<>
    inline const LetterValue CODON<33, 0u, 2u, 0u> = 15u;

    template<>
    inline const LetterValue CODON<33, 0u, 2u, 2u> = 8u;

    template<>
    inline const LetterValue CODON<33, 3u, 0u, 0u> = 19u;

    template<>
    inline const LetterValue CODON<33, 3u, 2u, 0u> = 18u;
}
