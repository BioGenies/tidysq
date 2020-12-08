#pragma once 

#include "tidysq/ops/remove_on_condition.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
        class OperationRemoveNA : public OperationRemoveOnCondition<INTERNAL_IN, INTERNAL_OUT> {
        public:
            OperationRemoveNA(const Alphabet &alphabet,
                              const bool by_letter) :
                    OperationRemoveOnCondition<INTERNAL_IN, INTERNAL_OUT>(
                            alphabet,
                            alphabet,
                            [&](const LetterValue &value) {
                                return this->OperationRemoveOnCondition<INTERNAL_IN, INTERNAL_OUT>::alph_.NA_value() != value;
                            },
                            by_letter) {};
        };
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sq<INTERNAL_OUT> remove_NA(const Sq<INTERNAL_IN> &sq,
                               const bool by_letter) {
        return sqapply(sq, ops::OperationRemoveNA<INTERNAL_IN, INTERNAL_OUT>(sq.alphabet(), by_letter));
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sequence<INTERNAL_OUT> remove_NA(const Sequence<INTERNAL_IN> &sequence,
                                     const Alphabet &alphabet,
                                     const bool by_letter) {
        return ops::OperationRemoveOnCondition<INTERNAL_IN, INTERNAL_OUT>(alphabet, by_letter)(sequence);
    }
}