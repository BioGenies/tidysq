#pragma once 

#include "tidysq/ops/Operation.h"
#include "tidysq/sqapply.h"
#include "tidysq/constants/complement_tables.h"

namespace tidysq {
    namespace internal {
        constexpr LetterValue read_complement(const SqType &type, const LetterValue& value) {
            switch (type) {
                case DNA_BSC:
                case RNA_BSC:
                    switch (value) {
                        case 0u:    return constants::COMPLEMENT<true, 0u>;
                        case 1u:    return constants::COMPLEMENT<true, 1u>;
                        case 2u:    return constants::COMPLEMENT<true, 2u>;
                        case 3u:    return constants::COMPLEMENT<true, 3u>;
                        case 4u:    return constants::COMPLEMENT<true, 4u>;
                        default:    return constants::COMPLEMENT<true, 7u>;
                    }
                case DNA_EXT:
                case RNA_EXT:
                    switch (value) {
                        case 0u:    return constants::COMPLEMENT<false, 0u>;
                        case 1u:    return constants::COMPLEMENT<false, 1u>;
                        case 2u:    return constants::COMPLEMENT<false, 2u>;
                        case 3u:    return constants::COMPLEMENT<false, 3u>;
                        case 4u:    return constants::COMPLEMENT<false, 4u>;
                        case 5u:    return constants::COMPLEMENT<false, 5u>;
                        case 6u:    return constants::COMPLEMENT<false, 6u>;
                        case 7u:    return constants::COMPLEMENT<false, 7u>;
                        case 8u:    return constants::COMPLEMENT<false, 8u>;
                        case 9u:    return constants::COMPLEMENT<false, 9u>;
                        case 10u:    return constants::COMPLEMENT<false, 10u>;
                        case 11u:    return constants::COMPLEMENT<false, 11u>;
                        case 12u:    return constants::COMPLEMENT<false, 12u>;
                        case 13u:    return constants::COMPLEMENT<false, 13u>;
                        case 14u:    return constants::COMPLEMENT<false, 14u>;
                        case 15u:    return constants::COMPLEMENT<false, 15u>;
                        default:    return constants::COMPLEMENT<false, 31u>;
                    }
                default:
                    throw std::invalid_argument("complement makes sense only for DNA and RNA sequences");
            }
        }
    }

    namespace ops {
        template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
        class OperationComplement : public OperationSqToSq<INTERNAL_IN, INTERNAL_OUT> {
            const AlphSize alph_size_;
            const SqType& type_;

        public:
            explicit OperationComplement(const AlphSize alph_size, const SqType &type) :
                    alph_size_(alph_size),
                    type_(type) {};

            void operator()(const Sequence<INTERNAL_IN> &sequence_in, Sequence<INTERNAL_OUT> &sequence_out) override {
                auto in_seq_iter = sequence_in.cbegin(alph_size_);
                auto out_seq_iter = sequence_out.begin(alph_size_);
                for (;  out_seq_iter != sequence_out.end(alph_size_) || in_seq_iter != sequence_in.cend(alph_size_);
                        ++in_seq_iter, ++out_seq_iter) {
                    out_seq_iter.assign(internal::read_complement(type_, *in_seq_iter));
                }
            }

            inline Sequence<INTERNAL_OUT> operator() (const Sequence<INTERNAL_IN> &sequence_in) override {
                //TODO: issue #57
                Sequence<INTERNAL_OUT> sequence_out = OperationSqToSq<INTERNAL_IN, INTERNAL_OUT>::initialize_element_out(sequence_in);
                operator()(sequence_in, sequence_out);
                return sequence_out;
            }
        };
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sq<INTERNAL_OUT> complement(const Sq<INTERNAL_IN> &sq) {
        return sqapply(sq, ops::OperationComplement<INTERNAL_IN, INTERNAL_OUT>(sq.alphabet().alphabet_size(), sq.type()));
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sequence<INTERNAL_OUT> complement(const Sequence<INTERNAL_IN> &sequence, const AlphSize alph_size, const SqType &type) {
        return ops::OperationComplement<INTERNAL_IN, INTERNAL_OUT>(alph_size, type)(sequence);
    }
}