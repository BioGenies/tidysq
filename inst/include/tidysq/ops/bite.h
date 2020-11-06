#pragma once

#include "tidysq/types/Sq.h"

namespace tidysq {
    template<InternalType INTERNAL>
    Sequence<INTERNAL> bite(const Sequence<INTERNAL> &sequence,
                      const std::vector<int> &indices,
                      const AlphSize &alph_size,
                      bool* warning_called) {
        Sequence<INTERNAL> out_sequence(
                internal::calculate_packed_internal_length(indices.size(), alph_size),
                indices.size()
        );
        const ElementPacked NA_value = 0xffu >> (8u - alph_size);

        auto index_iter = indices.begin();
        auto sequence_iter = sequence.begin(alph_size);
        auto out_sequence_iter = out_sequence.begin(alph_size);

        while (index_iter != indices.end() || out_sequence_iter != out_sequence.end(alph_size)) {
            ElementPacked element = NA_value;
            if (*index_iter <= sequence.originalLength()) {
                element = sequence_iter.access(*index_iter - 1);
            } else {
                *warning_called = true;
            }
            out_sequence_iter.assign(element);

            ++index_iter;
            ++out_sequence_iter;
        }
        return out_sequence;
    }

    template<InternalType INTERNAL>
    Sq<INTERNAL> bite(const SequenceIterator<INTERNAL> &it, const std::vector<int> &indices) {
        return bite(it.sequence_, indices, it.alph_size_, false);
    }

    template<InternalType INTERNAL>
    std::pair<std::string, Sq<INTERNAL>> bite(const Sq<INTERNAL> &sq, const std::vector<int> &indices) {
        Sq<INTERNAL> ret(sq.length(), sq.alphabet());
        bool warning_called = false;
        // TODO: replace with NULL once it works
        std::string NA_warning;

        for (LenSq i = 0; i < sq.length(); ++i) {
            // TODO: repair; why get()?
            ret[i] = bite(sq[i].get(), indices, sq.alphabet().alphabet_size(), &warning_called);
        }
        if (warning_called)
            NA_warning = "some sequences are subsetted with index bigger than length - NA introduced";

        return std::make_pair(NA_warning, ret);
    }
}
