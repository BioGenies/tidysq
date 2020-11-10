#pragma once

#include "tidysq/exports.h"

namespace tidysq {
    template<InternalType INTERNAL>
    std::vector<std::vector<Letter>> find_invalid_letters(const Sq<INTERNAL> &sq,
                                                          const SqType &type) {
        const Alphabet &alph = sq.alphabet();
        const Alphabet dest_alph = Alphabet(type, alph.NA_letter());

        std::vector<LetterValue> invalid_indices;
        for (LetterValue i = 0; i < alph.length(); ++i) {
            // TODO: would be nice to unfriend Alphabet and simply have an AlphabetIterator
            if (std::none_of(dest_alph.letters_.cbegin(), dest_alph.letters_.cend(),
                             [alph, i](const Letter& letter){ return alph[i] == letter; })) {
                invalid_indices.push_back(i);
            }
        }

        std::vector<std::vector<Letter>> ret;

        for (LenSq i = 0; i < sq.length(); ++i) {
            const Sequence<RCPP> sequence = sq[i];
            std::vector<Letter> invalid_found;
            for (const LetterValue &index : invalid_indices) {
                if (std::any_of(sequence.cbegin(alph.alphabet_size()), sequence.cend(alph.alphabet_size()),
                                [=](const ElementPacked elem){ return elem == index; })) {
                    invalid_found.push_back(alph[index]);
                }
            }
            ret.push_back(invalid_found);
        }

        return ret;
    }
}
