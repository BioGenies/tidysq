#pragma once 

#include "tidysq/ProtoSequence.h"
#include "tidysq/ops/Operation.h"
#include "tidysq/ops/pack.h"
#include "tidysq/util/random.h"
#include "tidysq/sqapply.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_OUT>
        class OperationRandomSq :
                public OperationVectorToVector<std::vector<LenSq>, LenSq,
                        Sq<INTERNAL_OUT>, Sequence<INTERNAL_OUT>> {
            const Alphabet &alphabet_;
            const bool &use_gap_;
            const std::vector<LetterValue> letter_values_;

            std::vector<LetterValue> filter_valid_letter_values() {
                std::vector<LetterValue> letter_values;
                for (const auto& pair : alphabet_) {
                    if (!((alphabet_.type() == AMI_BSC || alphabet_.type() == AMI_EXT) && pair.second == "*") &&
                        (use_gap_ || pair.second != "-")) {
                        letter_values.push_back(pair.first);
                    }
                }
                return letter_values;
            }

        public:
            OperationRandomSq(const Alphabet &alphabet,
                              const bool &use_gap) :
                    alphabet_(alphabet),
                    use_gap_(use_gap),
                    letter_values_(filter_valid_letter_values()) {};

            inline Sq<INTERNAL_OUT> initialize_vector_out(const std::vector<LenSq> &lengths, const LenSq from, const LenSq to) override {
                return Sq<INTERNAL_OUT>(to - from, alphabet_);
            }

            inline Sequence<INTERNAL_OUT> initialize_element_out(const LenSq &length) override {
                return util::reserve_space_for_packed<INTERNAL_OUT>(length, alphabet_.alphabet_size());
            }

            inline void operator() (const LenSq &length,
                                    Sequence<INTERNAL_OUT> &sequence) override {
                for (auto it = sequence.begin(alphabet_.alphabet_size()); it != sequence.end(alphabet_.alphabet_size()); ++it) {
                    it.assign(letter_values_[util::random_value<INTERNAL_OUT>(letter_values_.size())]);
                }
            }

            inline Sequence<INTERNAL_OUT> operator() (const LenSq &length) override {
                Sequence<INTERNAL_OUT> sequence = initialize_element_out(length);
                operator() (length, sequence);
                return sequence;
            }
        };
    }

    template<typename INTERNAL>
    Sq<INTERNAL> random_sq(const std::vector<LenSq> &lengths, const Alphabet &alphabet, const bool &use_gap) {
        return sqapply(lengths, ops::OperationRandomSq<INTERNAL>(alphabet, use_gap));
    }

    template<typename INTERNAL>
    Sq<INTERNAL> random_sq(const LenSq &n, const LenSq &len, const Alphabet &alphabet, const bool &use_gap) {
        std::vector<LenSq> lengths(n);
        std::fill(lengths.begin(), lengths.end(), len);
        return random_sq<INTERNAL>(lengths, alphabet, use_gap);
    }

    template<typename INTERNAL>
    Sq<INTERNAL> random_sequence(const LenSq length, const Alphabet &alphabet, const bool &use_gap) {
        return ops::OperationRandomSq<INTERNAL>(alphabet, use_gap)(length);
    }
}