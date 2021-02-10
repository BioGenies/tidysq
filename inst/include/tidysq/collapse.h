#pragma once

#include "tidysq/Sq.h"

namespace tidysq {
    template<typename INTERNAL>
    Sq<INTERNAL> collapse(const Sq<INTERNAL> &sq) {
        const Alphabet& alph = sq.alphabet();
        const AlphSize& alph_size = alph.alphabet_size();
        Sq<INTERNAL> ret(1, alph);

        // Early return in case there is exactly one sequence, so that no processing is done
        if (sq.size() == 1) {
            ret[0] = sq[0];
            return ret;
        }

        // First total sequence length (number of elements) is computed
        LenSq element_count = 0;
        for (LenSq i = 0; i < sq.size(); ++i) {
            element_count += sq[i].get().original_length();
        }
        Sequence<INTERNAL> sequence_out = Sequence<INTERNAL>((alph_size * element_count + 7) / 8, element_count);

        // Next an iterator is created and sequences are inputted one by one
        auto out_seq_iter = sequence_out.begin(alph_size);
        for (LenSq i = 0; i < sq.size(); ++i) {
            const Sequence<INTERNAL>& sequence_in = sq[i];
            auto in_seq_iter = sequence_in.cbegin(alph_size);
            while (out_seq_iter != sequence_out.end(alph_size) && in_seq_iter != sequence_in.cend(alph_size)) {
                out_seq_iter.assign(*in_seq_iter);
                ++in_seq_iter;
                ++out_seq_iter;
            }
        }
        ret[0] = sequence_out;

        return ret;
    }
}
