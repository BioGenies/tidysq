#pragma once

#include "tidysq/Sq.h"

namespace tidysq {
    typedef std::unordered_map<LetterValue, const LetterValue> ComplementTable;

    const ComplementTable bsc_complement_table = {
            {0u, 3u}, {1u, 2u}, {2u, 1u}, {3u, 0u}
    };

    const ComplementTable ext_complement_table = {
            {0u, 3u}, {1u, 2u}, {2u, 1u}, {3u, 0u},
            {4u, 4u}, {5u, 5u}, {6u, 7u}, {7u, 6u}, {8u, 9u}, {9u, 8u},
            {10u, 13u}, {11u, 12u}, {12u, 11u}, {13u, 10u},
            {14u, 14u}
    };

    template<typename INTERNAL>
    Sequence<INTERNAL> complement(const Sequence<INTERNAL> &sequence,
                                  const AlphSize &alph_size,
                                  const ComplementTable &table) {
        Sequence<INTERNAL> ret =
                util::reserve_space_for_packed<INTERNAL>(sequence.original_length(), alph_size);

        auto in_seq_iter = sequence.cbegin(alph_size);
        auto out_seq_iter = ret.begin(alph_size);
        while (out_seq_iter != ret.end(alph_size) || in_seq_iter != sequence.cend(alph_size)) {
            LetterValue in_letter = *in_seq_iter;
            if (table.count(in_letter) > 0) {
                out_seq_iter.assign(table.at(in_letter));
            } else {
                out_seq_iter.assign(in_letter);
            }
            ++in_seq_iter;
            ++out_seq_iter;
        }
        return ret;
    }

    template<typename INTERNAL>
    Sq<INTERNAL> complement(const Sq<INTERNAL> &sq) {
        Alphabet alph = sq.alphabet();
        Sq<INTERNAL> ret(sq.length(), alph);
        ComplementTable table;
        switch (sq.type()) {
            case DNA_BSC:
            case RNA_BSC:
                table = bsc_complement_table;
                break;
            case DNA_EXT:
            case RNA_EXT:
                table = ext_complement_table;
                break;
            default:
                throw std::invalid_argument("complement makes sense only for DNA and RNA sequences");
        }

        for (LenSq i = 0; i < sq.length(); ++i) {
            ret[i] = complement<INTERNAL>(sq[i].get(), alph.alphabet_size(), table);
        }
        return ret;
    }
}
