#pragma once

#include "tidysq/Sq.h"

namespace tidysq {
    template<typename INTERNAL>
    Sequence<INTERNAL> random_sequence(const int &len,
                                       const std::vector<LetterValue> &letter_values,
                                       const AlphSize &alph_size) {
        Sequence<INTERNAL> ret(util::calculate_packed_internal_length(len, alph_size), len);
        for (auto it = ret.begin(alph_size); it != ret.end(alph_size); ++it) {
            it.assign(letter_values[rand() % letter_values.size()]);
        }
        return ret;
    }

    template<typename INTERNAL>
    Sq<INTERNAL> random_sq(const std::vector<int> &len, const Alphabet &alphabet, const bool &use_gap) {
        Sq<INTERNAL> ret(len.size(), alphabet);
        std::vector<LetterValue> letter_values;
        for (std::pair<LetterValue, Letter> pair : alphabet) {
            if (!((alphabet.type() == AMI_BSC || alphabet.type() == AMI_EXT) && pair.second == "*") &&
                    (use_gap || pair.second != "-")) {
                letter_values.push_back(pair.first);
            }
        }

        for (LenSq i = 0; i < ret.length(); ++i) {
            ret[i] = random_sequence<INTERNAL>(len[i], letter_values, alphabet.alphabet_size());
        }
        return ret;
    }

    template<typename INTERNAL>
    Sq<INTERNAL> random_sq(const int &n, const int &len, const Alphabet &alphabet, const bool &use_gap) {
        std::vector<int> len_vector(n);
        std::fill(len_vector.begin(), len_vector.end(), len);
        return random_sq<INTERNAL>(len_vector, alphabet, use_gap);
    }
}
