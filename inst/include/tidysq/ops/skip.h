#pragma once

#include <tidysq/util/calculate_length.h>
#include "tidysq/ops/Operation.h"
#include "tidysq/sqapply.h"
#include "tidysq/util/ResultWrapper.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
        class OperationSkip : public OperationSqToSq<INTERNAL_IN, INTERNAL_OUT> {
            const AlphSize alphabet_size_;
            const std::vector<LenSq> indices_;

            [[nodiscard]] inline LenSq calculate_removed_indices_count(const Sequence<INTERNAL_IN> &sequence) const {
                // We have to count how many of these indices are actually influential
                return std::count_if(indices_.cbegin(), indices_.cend(), [=](auto &index) {
                    return index < sequence.original_length();
                });
            }

            [[nodiscard]] inline std::vector<LenSq> select_unique_indices(const std::vector<LenSq> &indices) {
                // First extract unique indices (as multiple instances of the same index doesn't change the result)
                std::vector<LenSq> unique_indices = indices;
                std::sort(unique_indices.begin(), unique_indices.end());
                auto last = std::unique(unique_indices.begin(), unique_indices.end());
                unique_indices.erase(last, unique_indices.end());
                return unique_indices;
            }

        public:
            explicit OperationSkip(const AlphSize alphabet_size_,
                                   const std::vector<LenSq> &indices) :
                    alphabet_size_(alphabet_size_),
                    indices_(select_unique_indices(indices)) {};


            inline Sequence<INTERNAL_OUT> initialize_element_out(const Sequence<INTERNAL_IN> &sequence) override {
                return util::reserve_space_for_packed<INTERNAL_OUT>(
                        sequence.original_length() - calculate_removed_indices_count(sequence),
                        alphabet_size_);
            }

            void operator()(const Sequence<INTERNAL_IN> &sequence_in, Sequence<INTERNAL_OUT> &sequence_out) override {
                auto in_sequence_iter = sequence_in.begin(alphabet_size_);
                auto out_sequence_iter = sequence_out.begin(alphabet_size_);

                while (out_sequence_iter != sequence_out.end(alphabet_size_)) {
                    if (std::none_of(indices_.cbegin(), indices_.cend(), [&in_sequence_iter](const long long int &index) {
                        return in_sequence_iter.index() == index;
                    })) {
                        out_sequence_iter.assign(*in_sequence_iter);
                        ++out_sequence_iter;
                    }
                    ++in_sequence_iter;
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
    Sq<INTERNAL_OUT> skip(const Sq<INTERNAL_IN> &sq, const std::vector<LenSq> &indices) {
        return sqapply(sq, ops::OperationSkip<INTERNAL_IN, INTERNAL_OUT>(sq.alphabet().alphabet_size(), indices));
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sequence<INTERNAL_OUT> skip(const Sequence<INTERNAL_IN> &sequence,
                                const AlphSize alphabet_size,
                                const std::vector<LenSq> &indices) {
        return ops::OperationSkip<INTERNAL_IN, INTERNAL_OUT>(alphabet_size, indices)(sequence);
    }
}