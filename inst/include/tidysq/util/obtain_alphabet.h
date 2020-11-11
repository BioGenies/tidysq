#pragma once

#include "tidysq/internal/obtain_alphabet.h"

namespace tidysq::util {
    template<InternalType INTERNAL>
    Alphabet obtain_alphabet(const typename TypeMapper<INTERNAL, STRING>::ProtoSqContentType &x,
                             const LenSq sample_size,
                             const Letter &NA_letter) {

        std::set<Letter> letters;

        if (NA_letter.length() == 0) {
            throw std::invalid_argument("'NA_letter' should have at least one character!");
        } else if (NA_letter.length() == 1) {
            letters = tidysq::internal::obtain_alphabet<INTERNAL, true>(x, sample_size, NA_letter);
        } else {
            letters = tidysq::internal::obtain_alphabet<INTERNAL, false>(x, sample_size, NA_letter);
        }

        return Alphabet(convert_set_to_vector(letters),UNT, NA_letter);
    }
}