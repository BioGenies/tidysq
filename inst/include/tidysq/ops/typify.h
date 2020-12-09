#pragma once

#include <utility>

#include "tidysq/ops/Operation.h"
#include "tidysq/sqapply.h"
#include "tidysq/ops/unpack.h"
#include "tidysq/ops/pack.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
        class OperationTypify : public OperationSqToSq<INTERNAL_IN, INTERNAL_OUT> {
        protected:
            const Alphabet alph_;
            const Alphabet dest_alph_;

        public:
            OperationTypify(Alphabet alphabet,
                            const SqType &dest_type) :
                    alph_(std::move(alphabet)),
                    dest_alph_(Alphabet(dest_type)) {
                if (!std::all_of(alph_.cbegin(), alph_.cend(), [=](const std::pair<LetterValue, Letter> &entry) {
                    return dest_alph_.contains(entry.second);
                })) {
                    throw std::invalid_argument("sq object contains letters that do not appear in the alphabet of target type");
                }
            };

            [[nodiscard]] Alphabet map_alphabet(const Alphabet &alphabet_in) const override {
                return dest_alph_;
            }

            bool may_return_early(const Sq<INTERNAL_IN> &vector_in) override {
                return std::is_same_v<INTERNAL_IN, INTERNAL_OUT> && alph_ == dest_alph_;
            }

            Sq<INTERNAL_OUT> return_early(const Sq<INTERNAL_IN> &vector_in) override {
                return  Sq<INTERNAL_OUT>(vector_in);
            }

            //this cannot be easily parallelized, as we do not know in advance the size of the resulting vector
            Sequence<INTERNAL_OUT> initialize_element_out(const Sequence<INTERNAL_IN> &sequence_in) override {
                return sequence_in;
            }

            void operator()(const Sequence<INTERNAL_IN> &sequence_in, Sequence<INTERNAL_OUT> &sequence_out) override {
                sequence_out = pack<STD_IT, STRING_PT, INTERNAL_OUT>(unpack<INTERNAL_IN, STD_IT, STRING_PT>(sequence_in, alph_), dest_alph_);
            }

            inline Sequence<INTERNAL_OUT> operator() (const Sequence<INTERNAL_IN> &sequence_in) override {
                Sequence<INTERNAL_OUT> sequence_out = initialize_element_out(sequence_in);
                operator()(sequence_in, sequence_out);
                return sequence_out;
            }
        };
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sq<INTERNAL_OUT> typify(const Sq<INTERNAL_IN> &sq,
                            const SqType &dest_type) {
        return sqapply(sq, ops::OperationTypify<INTERNAL_IN, INTERNAL_OUT>(sq.alphabet(), dest_type));
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sequence<INTERNAL_OUT> typify(const Sequence<INTERNAL_IN> &sequence,
                                  const Alphabet &alphabet,
                                  const SqType &dest_type) {
        return ops::OperationTypify<INTERNAL_IN, INTERNAL_OUT>(alphabet, dest_type)(sequence);
    }
}