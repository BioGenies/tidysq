#pragma once

#include "tidysq/internal/obtain_alphabet.h"

namespace tidysq::util {
    template<typename INTERNAL>
    Alphabet obtain_alphabet(const typename TypeBinder<INTERNAL, STRING_PT>::ProtoSqListConstructorType &x,
                             const LenSq sample_size,
                             const Letter &NA_letter,
                             const bool ignore_case) {

        std::set<Letter> letters;

        if (NA_letter.size() == 0) {
            throw std::invalid_argument("'NA_letter' should have at least one character!");
        } else if (NA_letter.size() == 1) {
            letters = tidysq::internal::obtain_alphabet<INTERNAL, true>(x, sample_size, NA_letter, ignore_case);
        } else {
            letters = tidysq::internal::obtain_alphabet<INTERNAL, false>(x, sample_size, NA_letter, ignore_case);
        }

        return Alphabet(convert_set_to_vector(letters),UNT, NA_letter);
    }
}