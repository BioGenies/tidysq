#pragma once

#include "tidysq/ops/Operation.h"
#include "tidysq/sqapply.h"
#include "tidysq/ProtoSequence.h"
#include "tidysq/ops/unpack.h"
#include "tidysq/ops/pack.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
        class OperationSubstituteLetters : public OperationSqToSq<INTERNAL_IN, INTERNAL_OUT> {
            const Alphabet alphabet_;
            const std::unordered_map<Letter, Letter> &encoding_;
            const Alphabet dest_alphabet_;
            const bool need_repacking_;

            Alphabet check_dest_alphabet() {
                // TODO: move the chunk below somewhere into Alphabet class
                std::vector<Letter> dest_letters;
                // We cannot rely on range-based loop, because it isn't ordered alphabetical by key. Thus it's necessary
                // to do manual looping. It should work as long as Alphabet class doesn't change to much.
                // TODO: come up with an idea for ordered looping that doesn't rely on this exact implementation of Alphabet
                for (LetterValue key = 0; key < alphabet_.size(); ++key) {
                    const Letter letter = alphabet_[key];
                    if (encoding_.count(letter) == 0) {
                        if (std::none_of(dest_letters.begin(), dest_letters.end(), [&](const Letter &other) {
                            return letter == other;
                        })) {
                            dest_letters.push_back(letter);
                        }
                    } else if (std::none_of(dest_letters.begin(), dest_letters.end(), [&](const Letter &other) {
                        return encoding_.at(letter) == other;
                    })) {
                        dest_letters.push_back(encoding_.at(letter));
                    }
                }
                return Alphabet(dest_letters, ATP, alphabet_.NA_letter());
            }

            bool check_need_repacking() {
                std::vector<Letter> values;
                std::set<Letter> unique_values;
                for (const auto& entry : encoding_) {
                    values.push_back(entry.second);
                    unique_values.insert(entry.second);
                }

                return std::any_of(values.cbegin(), values.cend(), [=](const Letter &value) {
                    return alphabet_.contains(value);
                }) || values.size() != unique_values.size();
            }

        public:
            explicit OperationSubstituteLetters(const Alphabet &alphabet,
                                                const std::unordered_map<Letter, Letter> &encoding) :
                    alphabet_(alphabet),
                    encoding_(encoding),
                    dest_alphabet_(check_dest_alphabet()),
                    need_repacking_(check_need_repacking()) {};


            [[nodiscard]] Alphabet map_alphabet(const Alphabet &alphabet_in) const override {
                return dest_alphabet_;
            }

            void operator()(const Sequence<INTERNAL_IN> &sequence_in, Sequence<INTERNAL_OUT> &sequence_out) override {
                if (need_repacking_) {
                    ProtoSequence<STD_IT, STRINGS_PT> unpacked = unpack<INTERNAL_IN, STD_IT, STRINGS_PT>(sequence_in, alphabet_);

                    // We have content as vector of strings, so that it's easier to swap them with encoding map
                    ProtoSequence<STD_IT, STRINGS_PT> encoded(unpacked.content());
                    for (LenSq index = 0; index < encoded.size(); ++index) {
                        Letter letter = alphabet_[sequence_in[{index, alphabet_.alphabet_size()}]];
                        if (encoding_.count(letter) > 0) {
                            encoded[index] = encoding_.at(letter);
                        }
                    }

                    (ops::OperationPack<STD_IT, STRINGS_PT, INTERNAL_OUT>(dest_alphabet_))(encoded, sequence_out);
                } else {
                    sequence_out = sequence_in;
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
    Sq<INTERNAL_OUT> substitute_letters(const Sq<INTERNAL_IN> &sq,
                                        const std::unordered_map<Letter, Letter> &encoding) {
        return sqapply(sq, ops::OperationSubstituteLetters<INTERNAL_IN, INTERNAL_OUT>(sq.alphabet(), encoding));
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sequence<INTERNAL_OUT> substitute_letters(const Sequence<INTERNAL_IN> &sequence,
                                              const Alphabet &alphabet,
                                              const std::unordered_map<Letter, Letter> &encoding) {
        return ops::OperationSubstituteLetters<INTERNAL_IN, INTERNAL_OUT>(alphabet, encoding)(sequence);
    }
}