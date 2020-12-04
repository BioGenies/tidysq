#pragma once

#include "tidysq/util/calculate_length.h"
#include "tidysq/ops/Operation.h"
#include "tidysq/sqapply.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
        class OperationReverse : public OperationSqToSq<INTERNAL_IN, INTERNAL_OUT> {
            const AlphSize alph_size_;
        public:
            explicit OperationReverse(AlphSize alph_size) :
                    alph_size_(alph_size) {};

            void operator()(const Sequence<INTERNAL_IN> &sequence_in, Sequence<INTERNAL_OUT> &sequence_out) override {
                // TODO: replace with const_reverse_iterator once implemented
                LenSq reverse_index = sequence_in.original_length() - 1;
                for (auto it = sequence_out.begin(alph_size_); it != sequence_out.end(alph_size_); ++it) {
                    it.assign(sequence_in[{reverse_index, alph_size_}]);
                    --reverse_index;
                }
            }
        };
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sq<INTERNAL_OUT> reverse(const Sq<INTERNAL_IN> &sq) {
        return sqapply(sq, ops::OperationReverse<INTERNAL_IN, INTERNAL_OUT>(sq.alphabet().alphabet_size()));
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sequence<INTERNAL_OUT> reverse(const Sequence<INTERNAL_IN> &sequence, const AlphSize alph_size) {
        return ops::OperationReverse<INTERNAL_IN, INTERNAL_OUT>(alph_size).
                template OperationVectorToVector<Sq<INTERNAL_IN>, Sequence<INTERNAL_IN>,
                Sq<INTERNAL_OUT>, Sequence<INTERNAL_OUT>>::operator()(sequence);
    }
}