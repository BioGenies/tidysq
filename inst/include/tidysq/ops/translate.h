#pragma once

#include <tidysq/util/calculate_length.h>
#include "tidysq/ops/Operation.h"
#include "tidysq/sqapply.h"

#include "tidysq/constants/translate_tables.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
        class OperationTranslate : public OperationSqToSq<INTERNAL_IN, INTERNAL_OUT> {
            const AlphSize NUC_BSC_ALPH_SIZE = Alphabet(DNA_BSC).alphabet_size();
            const AlphSize AMI_BSC_ALPH_SIZE = Alphabet(AMI_BSC).alphabet_size();

            const unsigned int table_;
            const bool interpret_as_stop_;

            inline unsigned int reduce_table(const unsigned int table) {
                // Some tables are actually identical to some others
                if (table == 7) return 4;
                if (table == 8) return 1;
                if (table == 11) return 1;
                return table;
            }

            // TODO: calculate codon tables at compile time - everything inside calculates with every codon!
            inline LetterValue codon_table(const LetterValue &codon_1,
                                           const LetterValue &codon_2,
                                           const LetterValue &codon_3) {
                // If this is non-standard table, then we have this kind of @Override
                if (table_ != 1) {
                    // The only way to include ambiguous translation tables (27, 28, 31)
                    // I don't really like this solution, maybe we should just drop the support for that
                    if ((table_ == 27 || table_ == 28 || table_ == 31) && !interpret_as_stop_) {
                        const auto& amb_codon_diff_table = constants::AMB_CODON_DIFF_TABLES.at(table_);
                        if (amb_codon_diff_table.count(codon_1) > 0 &&
                            amb_codon_diff_table.at(codon_1).count(codon_2) > 0 &&
                            amb_codon_diff_table.at(codon_1).at(codon_2).count(codon_3) > 0) {
                            return amb_codon_diff_table.at(codon_1).at(codon_2).at(codon_3);
                        }
                    }
                    // Then find correct table of differences
                    const auto& codon_diff_table = constants::CODON_DIFF_TABLES.at(table_);
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

        public:
            explicit OperationTranslate(const unsigned int &table,
                                        const bool &interpret_as_stop) :
                    table_(reduce_table(table)),
                    interpret_as_stop_(interpret_as_stop) {};

            [[nodiscard]] Alphabet map_alphabet(const Alphabet &alphabet_in) const override {
                return Alphabet(AMI_BSC, alphabet_in.NA_letter(), alphabet_in.ignores_case());
            }

            Sequence<INTERNAL_OUT> initialize_element_out(const Sequence<INTERNAL_IN> &sequence_in) override {
                return util::reserve_space_for_packed<INTERNAL_OUT>(sequence_in.original_length() / 3, AMI_BSC_ALPH_SIZE);
            }

            void operator()(const Sequence<INTERNAL_IN> &sequence_in, Sequence<INTERNAL_OUT> &sequence_out) override {
                if (sequence_out.size() > 0) {
                    auto input_it = sequence_in.cbegin(NUC_BSC_ALPH_SIZE);
                    auto output_it = sequence_out.begin(AMI_BSC_ALPH_SIZE);
                    while (input_it < sequence_in.cend(NUC_BSC_ALPH_SIZE) - 2) {
                        auto codon_1 = *input_it++;
                        auto codon_2 = *input_it++;
                        auto codon_3 = *input_it++;
                        output_it.assign(codon_table(codon_1, codon_2, codon_3));
                        ++output_it;
                    }
                }
            }

            inline Sequence<INTERNAL_OUT> operator() (const Sequence<INTERNAL_IN> &sequence_in) override {
                Sequence<INTERNAL_OUT> sequence_out = initialize_element_out(sequence_in);
                operator()(sequence_in, sequence_out);
                return sequence_out;
            }
        };
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sq<INTERNAL_OUT> translate(const Sq<INTERNAL_IN> &sq,
                               const unsigned int &table = 1,
                               const bool &interpret_as_stop = false) {
        return sqapply(sq, ops::OperationTranslate<INTERNAL_IN, INTERNAL_OUT>(table, interpret_as_stop));
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sequence<INTERNAL_OUT> translate(const Sequence<INTERNAL_IN> &sequence,
                                     const unsigned int &table = 1,
                                     const bool &interpret_as_stop = false) {
        return ops::OperationTranslate<INTERNAL_IN, INTERNAL_OUT>(table, interpret_as_stop)(sequence);
    }
}