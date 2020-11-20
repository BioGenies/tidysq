#pragma once

#include "tidysq/Sq.h"

namespace tidysq {
    template<typename INTERNAL>
    Sequence<INTERNAL> reverse(const Sequence<INTERNAL> &sequence,
                               const AlphSize &alph_size) {
        Sequence<INTERNAL> ret = util::reserve_space_for_packed<INTERNAL>(sequence.original_length(), alph_size);
        // TODO: replace with const_reverse_iterator once implemented
        LenSq reverse_index = sequence.original_length() - 1;
        for (auto it = ret.begin(alph_size); it != ret.end(alph_size); ++it) {
            it.assign(sequence[{reverse_index, alph_size}]);
            --reverse_index;
        }
        return ret;
    }

    template<typename INTERNAL>
    Sq<INTERNAL> reverse(const Sq<INTERNAL> &sq) {
        const Alphabet &alph = sq.alphabet();
        Sq<INTERNAL> ret(sq.length(), alph);

        for (LenSq i = 0; i < sq.length(); ++i) {
            ret[i] = reverse<INTERNAL>(sq[i].get(), alph.alphabet_size());
        }
        return ret;
    }
}
