#pragma once 

#include "tidysq/ops/Operation.h"
#include "tidysq/sqapply.h"
#include "tidysq/constants/complement_tables.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
        class OperationComplement : public OperationSqToSq<INTERNAL_IN, INTERNAL_OUT> {
            const AlphSize alph_size_;
            const internal::ComplementTable &table_;

            [[nodiscard]] const internal::ComplementTable & match_table(const SqType &type) const {
                switch (type) {
                    case DNA_BSC:
                    case RNA_BSC:
                        return constants::BSC_COMPLEMENT_TABLE;
                    case DNA_EXT:
                    case RNA_EXT:
                        return constants::EXT_COMPLEMENT_TABLE;
                    default:
                        throw std::invalid_argument("complement makes sense only for DNA and RNA sequences");
                }
            }

        public:
            explicit OperationComplement(const AlphSize alph_size, const SqType &type) :
                    alph_size_(alph_size),
                    table_(match_table(type)) {};

            void operator()(const Sequence<INTERNAL_IN> &sequence_in, Sequence<INTERNAL_OUT> &sequence_out) override {
                auto in_seq_iter = sequence_in.cbegin(alph_size_);
                auto out_seq_iter = sequence_out.begin(alph_size_);
                while (out_seq_iter != sequence_out.end(alph_size_) || in_seq_iter != sequence_in.cend(alph_size_)) {
                    LetterValue in_letter = *in_seq_iter;
                    if (table_.count(in_letter) > 0) {
                        out_seq_iter.assign(table_.at(in_letter));
                    } else {
                        out_seq_iter.assign(in_letter);
                    }
                    ++in_seq_iter;
                    ++out_seq_iter;
                }
            }
        };
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sq<INTERNAL_OUT> complement(const Sq<INTERNAL_IN> &sq) {
        return sqapply(sq, ops::OperationComplement<INTERNAL_IN, INTERNAL_OUT>(sq.alphabet().alphabet_size(), sq.type()));
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sequence<INTERNAL_OUT> complement(const Sequence<INTERNAL_IN> &sequence, const AlphSize alph_size, const SqType &type) {
        auto x = ops::OperationComplement<INTERNAL_IN, INTERNAL_OUT>(alph_size, type);
        return ops::OperationComplement<INTERNAL_IN, INTERNAL_OUT>(alph_size, type).
                template OperationVectorToVector<Sq<INTERNAL_IN>, Sequence<INTERNAL_IN>,
                        Sq<INTERNAL_OUT>, Sequence<INTERNAL_OUT>>::operator()(sequence);
    }
}