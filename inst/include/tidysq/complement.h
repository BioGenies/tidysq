#pragma once

#include "tidysq/Sq.h"
#include "tidysq/constants/complement_tables.h"

namespace tidysq {
    template<typename INTERNAL>
    Sequence<INTERNAL> complement(const Sequence<INTERNAL> &sequence,
                                  const AlphSize &alph_size,
                                  const internal::ComplementTable &table) {
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
        internal::ComplementTable table;
        switch (sq.type()) {
            case DNA_BSC:
            case RNA_BSC:
                table = constants::BSC_COMPLEMENT_TABLE;
                break;
            case DNA_EXT:
            case RNA_EXT:
                table = constants::EXT_COMPLEMENT_TABLE;
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
