#pragma once

#include "tidysq/ops/Operation.h"
#include "tidysq/sqapply.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN>
        class OperationFindInvalidLetters : public OperationVectorToVector<Sq<INTERNAL_IN>, Sequence<INTERNAL_IN>,
                                                        std::vector<std::vector<Letter>>, std::vector<Letter>> {
        private:
            const Alphabet alph_;
            const Alphabet dest_alph_;
            const std::vector<LetterValue> invalid_indices_;

            std::vector<LetterValue> find_invalid_indices() const {
                std::vector<LetterValue> invalid_indices = {};
                for (LetterValue i = 0; i < alph_.size(); ++i) {
                    if (std::none_of(dest_alph_.cbegin(), dest_alph_.cend(),
                                     [&](const auto& pair){ return alph_[i] == pair.second; })) {
                        invalid_indices.push_back(i);
                    }
                }
                return invalid_indices;
            }

        public:
            OperationFindInvalidLetters(Alphabet alphabet,
                                        const SqType &dest_type) :
                    alph_(std::move(alphabet)),
                    dest_alph_(Alphabet(dest_type)),
                    invalid_indices_(find_invalid_indices()) {};

            bool may_return_early(const Sq<INTERNAL_IN> &sq) override {
                return alph_ == dest_alph_;
            }

            std::vector<std::vector<Letter>> return_early(const Sq<INTERNAL_IN> &sq) override {
                return  std::vector<std::vector<Letter>>(sq.size());
            }

            std::vector<std::vector<Letter>> initialize_vector_out(
                    const Sq<INTERNAL_IN> &vector_in, const LenSq from, const LenSq to) override {
                return std::vector<std::vector<Letter>>(to - from);
            };


            std::vector<Letter> initialize_element_out(const Sequence<INTERNAL_IN> &sequence) override {
                return {};
            }

            void operator()(const Sequence<INTERNAL_IN> &sequence, std::vector<Letter> &letters) override {
                for (const LetterValue &index : invalid_indices_) {
                    if (std::any_of(sequence.cbegin(alph_.alphabet_size()), sequence.cend(alph_.alphabet_size()),
                                    [&](const ElementPacked elem){ return elem == index; })) {
                        letters.push_back(alph_[index]);
                    }
                }
            }

            inline std::vector<Letter> operator() (const Sequence<INTERNAL_IN> &sequence) override {
                std::vector<Letter> vec = initialize_element_out(sequence);
                operator()(sequence, vec);
                return vec;
            }
        };
    }

    template<typename INTERNAL_IN>
    std::vector<std::vector<Letter>> find_invalid_letters(const Sq<INTERNAL_IN> &sq,
                                                          const SqType &dest_type) {
        return sqapply(sq, ops::OperationFindInvalidLetters<INTERNAL_IN>(sq.alphabet(), dest_type));
    }

    template<typename INTERNAL_IN>
    std::vector<Letter> find_invalid_letters(const Sequence<INTERNAL_IN> &sequence,
                                             const Alphabet &alphabet,
                                             const SqType &dest_type) {
        return ops::OperationFindInvalidLetters<INTERNAL_IN>(alphabet, dest_type)(sequence);
    }
}