#pragma once

#include <map>

#include "tidysq/constants/translate_tables.h"

namespace tidysq {
    inline LetterValue codon_table(int table,
                                   const LetterValue &codon_1,
                                   const LetterValue &codon_2,
                                   const LetterValue &codon_3,
                                   const bool &interpret_as_stop) {
        // Some tables are actually identical to some others
        if (table == 7) table = 4;
        if (table == 8) table = 1;
        if (table == 11) table = 1;
        // If this is non-standard table, then we have this kind of @Override
        if (table != 1) {
            // The only way to include ambiguous translation tables (27, 28, 31)
            // I don't really like this solution, maybe we should just drop the support for that
            if ((table == 27 || table == 28 || table == 31) && !interpret_as_stop) {
                const auto& amb_codon_diff_table = constants::AMB_CODON_DIFF_TABLES.at(table);
                if (amb_codon_diff_table.count(codon_1) > 0 &&
                        amb_codon_diff_table.at(codon_1).count(codon_2) > 0 &&
                        amb_codon_diff_table.at(codon_1).at(codon_2).count(codon_3) > 0) {
                    return amb_codon_diff_table.at(codon_1).at(codon_2).at(codon_3);
                }
            }
            // Then find correct table of differences
            const auto& codon_diff_table = constants::CODON_DIFF_TABLES.at(table);
            if (codon_diff_table.count(codon_1) > 0 &&
                    codon_diff_table.at(codon_1).count(codon_2) > 0 &&
                    codon_diff_table.at(codon_1).at(codon_2).count(codon_3) > 0) {
                const auto amino_acid = codon_diff_table.at(codon_1).at(codon_2).at(codon_3);
                // TODO: handle case if (amino_acid == 31u) (i.e. NA_letter)
                return amino_acid;
            }
        }
        return constants::CODON_TABLE_1.at(codon_1).at(codon_2).at(codon_3);
    }

    template<typename INTERNAL>
    Sequence<INTERNAL> translate(const Sequence<INTERNAL> &sequence,
                                 const int &table,
                                 const bool &interpret_as_stop,
                                 const AlphSize &input_alph_size,
                                 const AlphSize &output_alph_size) {
        LenSq sequence_length = sequence.original_length() / 3;
        Sequence<INTERNAL> ret = util::reserve_space_for_packed<INTERNAL>(sequence_length, output_alph_size);

        if (sequence_length > 0) {
            auto input_it = sequence.cbegin(input_alph_size);
            auto output_it = ret.begin(output_alph_size);
            while (input_it < sequence.cend(input_alph_size) - 2) {
                auto codon_1 = *input_it++;
                auto codon_2 = *input_it++;
                auto codon_3 = *input_it++;
                output_it.assign(codon_table(table, codon_1, codon_2, codon_3, interpret_as_stop));
                ++output_it;
            }
        }
        return ret;
    }

    template<typename INTERNAL>
    Sq<INTERNAL> translate(const Sq<INTERNAL> &sq,
                           const int &table,
                           const bool &interpret_as_stop) {
        const Alphabet& alph = sq.alphabet();
        Sq<INTERNAL> ret(sq.length(), Alphabet(AMI_BSC, alph.NA_letter()));

        for (LenSq i = 0; i < sq.length(); ++i) {
            ret[i] = translate(sq[i].get(), table, interpret_as_stop,
                    alph.alphabet_size(), ret.alphabet().alphabet_size());
        }
        return ret;
    }
}
