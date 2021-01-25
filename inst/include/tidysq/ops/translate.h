#pragma once

#include <tidysq/util/calculate_length.h>
#include "tidysq/ops/Operation.h"
#include "tidysq/sqapply.h"

#include "tidysq/constants/translate_tables.h"

namespace tidysq {
    namespace internal {

#define THIRD_CODON(CODON_1, CODON_2)                                                                                  \
switch (value_3) {                                                                                                     \
    case 0u:    return constants::CODON<TABLE, CODON_1##u, CODON_2##u, 0u>;                                            \
    case 1u:    return constants::CODON<TABLE, CODON_1##u, CODON_2##u, 1u>;                                            \
    case 2u:    return constants::CODON<TABLE, CODON_1##u, CODON_2##u, 2u>;                                            \
    case 3u:    return constants::CODON<TABLE, CODON_1##u, CODON_2##u, 3u>;                                            \
    default:    throw std::invalid_argument("translation must be made with four standard DNA/RNA letters only");       \
}

#define SECOND_CODON(CODON_1)                                                                                          \
switch (value_2) {                                                                                                     \
    case 0u:    THIRD_CODON(CODON_1, 0)                                                                                \
    case 1u:    THIRD_CODON(CODON_1, 1)                                                                                \
    case 2u:    THIRD_CODON(CODON_1, 2)                                                                                \
    case 3u:    THIRD_CODON(CODON_1, 3)                                                                                \
    default:    throw std::invalid_argument("translation must be made with four standard DNA/RNA letters only");       \
}

#define FIRST_CODON                                                                                                    \
switch (value_1) {                                                                                                     \
    case 0u:    SECOND_CODON(0)                                                                                        \
    case 1u:    SECOND_CODON(1)                                                                                        \
    case 2u:    SECOND_CODON(2)                                                                                        \
    case 3u:    SECOND_CODON(3)                                                                                        \
    default:    throw std::invalid_argument("translation must be made with four standard DNA/RNA letters only");       \
}

        template<int TABLE>
        constexpr LetterValue read_codon(LetterValue value_1, LetterValue value_2, LetterValue value_3) {
            FIRST_CODON
        }

#undef FIRST_CODON
#undef SECOND_CODON
#undef THIRD_CODON

        constexpr LetterValue read_codon(int table, LetterValue value_1, LetterValue value_2, LetterValue value_3) {
            switch (table) {
                case 1:
                case 8:
                case 11:
                    return read_codon<1>(value_1, value_2, value_3);
                case 2:
                    return read_codon<2>(value_1, value_2, value_3);
                case 3:
                    return read_codon<3>(value_1, value_2, value_3);
                case 4:
                case 7:
                    return read_codon<4>(value_1, value_2, value_3);
                case 5:
                    return read_codon<5>(value_1, value_2, value_3);
                case 6:
                    return read_codon<6>(value_1, value_2, value_3);
                case 9:
                    return read_codon<9>(value_1, value_2, value_3);
                case 10:
                    return read_codon<10>(value_1, value_2, value_3);
                case 12:
                    return read_codon<12>(value_1, value_2, value_3);
                case 13:
                    return read_codon<13>(value_1, value_2, value_3);
                case 14:
                    return read_codon<14>(value_1, value_2, value_3);
                case 15:
                    return read_codon<15>(value_1, value_2, value_3);
                case 16:
                    return read_codon<16>(value_1, value_2, value_3);
                case 21:
                    return read_codon<21>(value_1, value_2, value_3);
                case 22:
                    return read_codon<22>(value_1, value_2, value_3);
                case 23:
                    return read_codon<23>(value_1, value_2, value_3);
                case 24:
                    return read_codon<24>(value_1, value_2, value_3);
                case 25:
                    return read_codon<25>(value_1, value_2, value_3);
                case 26:
                    return read_codon<26>(value_1, value_2, value_3);
                case 29:
                    return read_codon<29>(value_1, value_2, value_3);
                case 30:
                    return read_codon<30>(value_1, value_2, value_3);
                case 33:
                    return read_codon<33>(value_1, value_2, value_3);
                default:
                    throw std::invalid_argument("specified table doesn't exist");
            }
        }
    }

    namespace ops {
        template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
        class OperationTranslate : public OperationSqToSq<INTERNAL_IN, INTERNAL_OUT> {
            const AlphSize NUC_BSC_ALPH_SIZE = Alphabet(DNA_BSC).alphabet_size();
            const AlphSize AMI_BSC_ALPH_SIZE = Alphabet(AMI_BSC).alphabet_size();

            const unsigned int table_;

        public:
            explicit OperationTranslate(const unsigned int &table) :
                    table_(table) {};

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
                        output_it.assign(internal::read_codon(table_, codon_1, codon_2, codon_3));
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
                               const unsigned int &table = 1) {
        return sqapply(sq, ops::OperationTranslate<INTERNAL_IN, INTERNAL_OUT>(table));
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sequence<INTERNAL_OUT> translate(const Sequence<INTERNAL_IN> &sequence,
                                     const unsigned int &table = 1) {
        return ops::OperationTranslate<INTERNAL_IN, INTERNAL_OUT>(table)(sequence);
    }
}