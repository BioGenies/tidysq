#pragma once

#include <utility>

#include "tidysq/ops/Operation.h"
#include "tidysq/sqapply.h"
#include "tidysq/ops/unpack.h"
#include "tidysq/ops/pack.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
        class OperationRemoveOnCondition : public OperationSqToSq<INTERNAL_IN, INTERNAL_OUT> {
        protected:
            const Alphabet alph_;
            const Alphabet dest_alph_;
            const std::function<bool(const LetterValue &)> condition_;
            const bool by_letter_;

        public:
            OperationRemoveOnCondition(Alphabet alphabet,
                                       Alphabet dest_alphabet,
                                       std::function<bool(const LetterValue &)> condition,
                                       const bool by_letter) :
                    alph_(std::move(alphabet)),
                    dest_alph_(std::move(dest_alphabet)),
                    condition_(std::move(condition)),
                    by_letter_(by_letter) {};

            [[nodiscard]] Alphabet map_alphabet(const Alphabet &alphabet_in) const override {
                return dest_alph_;
            }

            //this cannot be easily parallelized, as we do not know in advance the size of the resulting vector
            Sequence<INTERNAL_OUT> initialize_element_out(const Sequence<INTERNAL_IN> &sequence_in) override {
                return Sequence<INTERNAL_OUT>{0,0};
            }

            void operator()(const Sequence<INTERNAL_IN> &sequence_in, Sequence<INTERNAL_OUT> &sequence_out) override {
                if (by_letter_) {
                    typename ProtoSequence<STD_IT, STRINGS_PT>::ContentStorageType selected_letters;
                    for (auto it = sequence_in.cbegin(alph_.alphabet_size()); it != sequence_in.cend(alph_.alphabet_size()); ++it) {
                        if (condition_(*it)) {
                            selected_letters.push_back(alph_[*it]);
                        }
                    }
                    sequence_out = pack<STD_IT, STRINGS_PT, INTERNAL_OUT>(ProtoSequence<STD_IT, STRINGS_PT>{selected_letters}, dest_alph_);
                } else {
                    if (std::all_of(sequence_in.cbegin(alph_.alphabet_size()), sequence_in.cend(alph_.alphabet_size()),
                                    [=](const LetterValue element) { return condition_(element); })) {
                        if (alph_ != dest_alph_) {
                            sequence_out = pack<STD_IT, STRINGS_PT, INTERNAL_OUT>(
                                    unpack<INTERNAL_IN, STD_IT, STRINGS_PT>(sequence_in, alph_),
                                    dest_alph_);
                        } else {
                            sequence_out = sequence_in;
                        }
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
    Sq<INTERNAL_OUT> remove_on_condition(const Sq<INTERNAL_IN> &sq,
                                         const Alphabet &dest_alphabet,
                                         const std::function<bool(const LetterValue &)> &condition,
                                         const bool by_letter) {
        return sqapply(sq, ops::OperationRemoveOnCondition<INTERNAL_IN, INTERNAL_OUT>(sq.alphabet(), dest_alphabet, condition, by_letter));
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sequence<INTERNAL_OUT> remove_on_condition(const Sequence<INTERNAL_IN> &sequence,
                                               const Alphabet &alphabet,
                                               const Alphabet &dest_alphabet,
                                               const std::function<bool(const LetterValue &)> &condition,
                                               const bool by_letter) {
        return ops::OperationRemoveOnCondition<INTERNAL_IN, INTERNAL_OUT>(alphabet, dest_alphabet, condition, by_letter)(sequence);
    }
}