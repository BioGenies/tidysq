#pragma once

#include "tidysq/Sq.h"
#include "tidysq/util/find_common_size.h"

namespace tidysq {
    template<typename INTERNAL>
    Sq<INTERNAL> paste(const std::vector<Sq<INTERNAL>> &list_of_sq) {
        // If there are no sq, then we don't know what SqType to return
        if (list_of_sq.empty()) {
            throw std::invalid_argument("cannot paste empty vector");
        }

        // Now we can assume that there exists an element under index 0
        const Alphabet& alph = list_of_sq[0].alphabet();
        const AlphSize& alph_size = alph.alphabet_size();

        if (std::any_of(list_of_sq.cbegin(), list_of_sq.cend(), [alph](const Sq<INTERNAL>& other) {
            return other.alphabet() != alph;
        })) {
            throw std::invalid_argument("pasting sq objects with different alphabets is not implemented");
        }

        const LenSq& common_size = util::find_common_size(list_of_sq);
        Sq<INTERNAL> ret(common_size, alph);

        // First total sequence lengths (number of elements) is computed
        std::vector<LenSq> element_counts(common_size);
        for (const Sq<INTERNAL>& sq : list_of_sq) {
            // If sq is a scalar (length 1), then increment all counts by its original_length
            if (sq.size() == 1) {
                for (LenSq& element_count : element_counts) {
                    element_count += sq[0].get().original_length();
                }
            } else {
                for (u_LenSq i = 0; i < element_counts.size(); ++i) {
                    element_counts[i] += sq[i].get().original_length();
                }
            }
        }

        for (LenSq i = 0; i < ret.size(); ++i) {
            Sequence<INTERNAL> sequence_out =
                    Sequence<INTERNAL>((alph_size * element_counts[i] + 7) / 8, element_counts[i]);
            auto out_seq_iter = sequence_out.begin(alph_size);
            for (const Sq<INTERNAL>& sq : list_of_sq) {
                // If sq is a scalar (length 1), then append the only Sequence present
                const Sequence<INTERNAL>& sequence_in = sq[sq.size() == 1 ? 0 : i];
                auto in_seq_iter = sequence_in.cbegin(alph_size);
                while (out_seq_iter != sequence_out.end(alph_size) && in_seq_iter != sequence_in.cend(alph_size)) {
                    out_seq_iter.assign(*in_seq_iter);
                    ++in_seq_iter;
                    ++out_seq_iter;
                }
            }
            ret[i] = sequence_out;
        }

        return ret;
    }
}
