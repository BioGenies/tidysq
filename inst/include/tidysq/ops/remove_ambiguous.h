#pragma once

#include "tidysq/ops/remove_on_condition.h"

namespace tidysq {
    namespace ops {
        template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
        class OperationRemoveAmbiguous : public OperationRemoveOnCondition<INTERNAL_IN, INTERNAL_OUT> {
            bool early_return_;

            Alphabet match_dest_alph(const Alphabet &alphabet) {
                switch (alphabet.type()) {
                    case AMI_BSC:
                        early_return_ = true;
                    case AMI_EXT:
                        return Alphabet(AMI_BSC);
                    case DNA_BSC:
                        early_return_ = true;
                    case DNA_EXT:
                        return Alphabet(DNA_BSC);
                    case RNA_BSC:
                        early_return_ = true;
                    case RNA_EXT:
                        return Alphabet(RNA_BSC);
                    default:
                        throw std::invalid_argument("sq object must have alphabet type that has corresponding extended alphabet type");
                }
            }

            bool may_return_early(const Sq<INTERNAL_IN> &vector_in) override {
                return early_return_;
            }



        public:
            OperationRemoveAmbiguous(const Alphabet &alphabet,
                                     const bool by_letter) :
                    OperationRemoveOnCondition<INTERNAL_IN, INTERNAL_OUT>(
                            alphabet,
                            match_dest_alph(alphabet),
                            [&](const LetterValue &value){
                                return this->OperationRemoveOnCondition<INTERNAL_IN, INTERNAL_OUT>::dest_alph_.contains(
                                this->OperationRemoveOnCondition<INTERNAL_IN, INTERNAL_OUT>::alph_[value]) ||
                                this->OperationRemoveOnCondition<INTERNAL_IN, INTERNAL_OUT>::alph_.NA_value() == value;
                            },
                            by_letter),
                    early_return_(false) {};
        };
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sq<INTERNAL_OUT> remove_ambiguous(const Sq<INTERNAL_IN> &sq,
                               const bool by_letter) {
        return sqapply(sq, ops::OperationRemoveAmbiguous<INTERNAL_IN, INTERNAL_OUT>(sq.alphabet(), by_letter));
    }

    template<typename INTERNAL_IN, typename INTERNAL_OUT = INTERNAL_IN>
    Sequence<INTERNAL_OUT> remove_ambiguous(const Sequence<INTERNAL_IN> &sequence,
                                            const Alphabet &alphabet,
                                            const bool by_letter) {
        return ops::OperationRemoveAmbiguous<INTERNAL_IN, INTERNAL_OUT>(alphabet, by_letter)(sequence);
    }
}