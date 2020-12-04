#pragma once

#include "tidysq/Sq.h"
#include "tidysq/util/calculate_length.h"

namespace tidysq {
    template<typename INTERNAL>
    Sequence<INTERNAL> bite(const Sequence<INTERNAL> &sequence,
                            const std::vector<long long int> &indices,
                            const AlphSize &alph_size,
                            bool* warning_called) {
        Sequence<INTERNAL> out_sequence(
                util::calculate_packed_internal_length(indices.size(), alph_size),
                indices.size()
        );
        const ElementPacked NA_value = 0xffu >> (8u - alph_size);

        auto index_iter = indices.begin();
        auto out_sequence_iter = out_sequence.begin(alph_size);

        while (out_sequence_iter != out_sequence.end(alph_size)) {
            ElementPacked element = NA_value;
            if (*index_iter < sequence.original_length()) {
                element = sequence[{*index_iter, alph_size}];
            } else {
                *warning_called = true;
            }
            out_sequence_iter.assign(element);

            ++index_iter;
            ++out_sequence_iter;
        }
        return out_sequence;
    }

    template<typename INTERNAL>
    Sequence<INTERNAL> bite(const typename Sequence<INTERNAL>::const_iterator &it,
                            const std::vector<long long int> &indices) {
        bool* warning_called = new bool;
        auto ret = bite(it.sequence_, indices, it.alph_size_, warning_called);
        delete warning_called;
        return ret;
    }

    template<typename INTERNAL>
    std::pair<std::string, Sq<INTERNAL>> bite(const Sq<INTERNAL> &sq,
                                              const std::vector<long long int> &indices) {
        Sq<INTERNAL> ret(sq.size(), sq.alphabet());
        bool warning_called = false;
        // TODO: replace with NULL once it works
        std::string NA_warning;

        for (LenSq i = 0; i < sq.size(); ++i) {
            ret[i] = bite(sq[i].get(), indices, sq.alphabet().alphabet_size(), &warning_called);
        }
        if (warning_called)
            NA_warning = "some sequences are subsetted with index bigger than length - NA introduced";

        return std::make_pair(NA_warning, ret);
    }

    template<typename INTERNAL>
    Sequence<INTERNAL> skip(const Sequence<INTERNAL> &sequence,
                            const std::vector<long long int> &indices,
                            const AlphSize &alph_size) {
        // We have to count how many of these indices are actually influential
        long long int removed_indices_count =
                std::count_if(indices.cbegin(), indices.cend(), [=](auto &index) {
                    return index < sequence.original_length();
                });
        Sequence<INTERNAL> out_sequence(
                util::calculate_packed_internal_length(sequence.original_length() - removed_indices_count, alph_size),
                sequence.original_length() - removed_indices_count
        );

        auto in_sequence_iter = sequence.begin(alph_size);
        auto out_sequence_iter = out_sequence.begin(alph_size);

        while (out_sequence_iter != out_sequence.end(alph_size)) {
            if (std::none_of(indices.cbegin(), indices.cend(), [&in_sequence_iter](const long long int &index) {
                return in_sequence_iter.index() == index;
            })) {
                out_sequence_iter.assign(*in_sequence_iter);
                ++out_sequence_iter;
            }
            ++in_sequence_iter;
        }
        return out_sequence;
    }

    template<typename INTERNAL>
    Sq<INTERNAL> skip(const Sq<INTERNAL> &sq,
                      const std::vector<long long int> &indices) {
        Sq<INTERNAL> ret(sq.size(), sq.alphabet());

        // First extract unique indices (as multiple instances of the same index doesn't change the result)
        std::vector<long long int> unique_indices = indices;
        std::sort(unique_indices.begin(), unique_indices.end());
        auto last = std::unique(unique_indices.begin(), unique_indices.end());
        unique_indices.erase(last, unique_indices.end());

        for (LenSq i = 0; i < sq.size(); ++i) {
            ret[i] = skip(sq[i].get(), unique_indices, sq.alphabet().alphabet_size());
        }
        return ret;
    }
}
