#pragma once

#include "tidysq/Sq.h"

namespace tidysq {
    template<typename INTERNAL>
    Sq<INTERNAL> typify(const Sq<INTERNAL> &sq,
                        const SqType &type) {
        // Early return whenever current type is equal to target type
        if (sq.type() == type) {
            return sq;
        }

        const Alphabet &alph = sq.alphabet();
        const Alphabet dest_alph = Alphabet(type, alph.NA_letter());

        // Input alphabet must be a subset of target alphabet, otherwise some letters cannot be encoded
        if (!std::all_of(alph.cbegin(), alph.cend(), [=](const std::pair<LetterValue, Letter> &entry) {
            return dest_alph.contains(entry.second);
        })) {
            throw std::invalid_argument("sq object contains letters that do not appear in the alphabet of target type");
        }

        ProtoSq<INTERNAL, STRINGS_PT> unpacked = sq.template unpack<INTERNAL, STRINGS_PT>();
        unpacked.alphabet() = dest_alph;
        return unpacked.template pack<INTERNAL>();
    }
}
