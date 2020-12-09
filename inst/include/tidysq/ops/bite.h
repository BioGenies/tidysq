#pragma once

#include <tidysq/util/calculate_length.h>
#include "tidysq/ops/Operation.h"
#include "tidysq/sqapply.h"
#include "tidysq/util/ResultWrapper.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
        class OperationBite : public OperationSqToSq<INTERNAL_IN, INTERNAL_OUT> {
            const AlphSize alphabet_size_;
            const ElementPacked NA_value_;
            const std::vector<LenSq> &indices_;
            bool warning_called_;

            [[nodiscard]] inline LetterValue calculate_NA_value() const {
                return 0xffu >> (8u - alphabet_size_);
            }

        public:
            explicit OperationBite(const AlphSize alphabet_size_,
                                   const std::vector<LenSq> &indices) :
                    alphabet_size_(alphabet_size_),
                    NA_value_(calculate_NA_value()),
                    indices_(indices),
                    warning_called_(false) {};


            inline Sequence<INTERNAL_OUT> initialize_element_out(const Sequence<INTERNAL_IN> &sequence) override {
                return util::reserve_space_for_packed<INTERNAL_OUT>(indices_.size(), alphabet_size_);
            }

            void operator()(const Sequence<INTERNAL_IN> &sequence_in, Sequence<INTERNAL_OUT> &sequence_out) override {
                auto index_iter = indices_.begin();
                auto out_sequence_iter = sequence_out.begin(alphabet_size_);

                while (out_sequence_iter != sequence_out.end(alphabet_size_)) {
                    ElementPacked element = NA_value_;
                    if (*index_iter < sequence_in.original_length()) {
                        element = sequence_in[{*index_iter, alphabet_size_}];
                    } else {
                        warning_called_ = true;
                    }
                    out_sequence_iter.assign(element);

                    ++index_iter;
                    ++out_sequence_iter;
                }
            }

            inline Sequence<INTERNAL_OUT> operator() (const Sequence<INTERNAL_IN> &sequence_in) override {
                Sequence<INTERNAL_OUT> sequence_out = initialize_element_out(sequence_in);
                operator()(sequence_in, sequence_out);
                return sequence_out;
            }

            std::optional<std::string> warning() {
                return warning_called_ ?
                std::make_optional("some sequences are subsetted with index bigger than length - NA introduced") :
                std::nullopt;
            }

        };
    }


    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    util::ResultWrapper<Sq<INTERNAL_OUT>> bite(const Sq<INTERNAL_IN> &sq, const std::vector<LenSq> &indices) {
        auto op = ops::OperationBite<INTERNAL_IN, INTERNAL_OUT>(sq.alphabet().alphabet_size(), indices);
        auto ret = sqapply(sq, op);
        return util::ResultWrapper(ret, op.warning());
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    util::ResultWrapper<Sequence<INTERNAL_OUT>> bite(const Sequence<INTERNAL_IN> &sequence,
                                                     const AlphSize alphabet_size,
                                                     const std::vector<LenSq> &indices) {
        auto op = ops::OperationBite<INTERNAL_IN, INTERNAL_OUT>(alphabet_size, indices);
        auto ret = op(sequence);
        return util::ResultWrapper(ret, op.warning());
    }
}