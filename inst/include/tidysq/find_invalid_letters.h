#pragma once

#include "tidysq/Alphabet.h"
#include "tidysq/Sequence.h"

namespace tidysq {
    template<typename INTERNAL>
    std::vector<std::vector<Letter>> find_invalid_letters(const Sq<INTERNAL> &sq,
                                                          const SqType &type) {
        const Alphabet &alph = sq.alphabet();
        const Alphabet dest_alph = Alphabet(type, alph.NA_letter());

        std::vector<LetterValue> invalid_indices;
        for (LetterValue i = 0; i < alph.size(); ++i) {
            if (std::none_of(dest_alph.cbegin(), dest_alph.cend(),
                             [alph, i](const auto& pair){ return alph[i] == pair.second; })) {
                invalid_indices.push_back(i);
            }
        }

        std::vector<std::vector<Letter>> ret;

        for (LenSq i = 0; i < sq.size(); ++i) {
            const Sequence<RCPP_IT> sequence = sq[i];
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
