#ifndef TIDYSQ_TYPES_STDSQPROTO_H
#define TIDYSQ_TYPES_STDSQPROTO_H

#include "types_StdSequence.h"
#include "impl_Operation.h"
#include "sqapply.h"

namespace tidysq {
    template<typename STD_SEQUENCE_PROTO>
    class StdSqProto {
        std::vector<STD_SEQUENCE_PROTO> content_;
        Alphabet alphabet_;
    public:
        typedef STD_SEQUENCE_PROTO SequenceType;
        typedef Alphabet AlphabetType;

        StdSqProto(std::vector<STD_SEQUENCE_PROTO> content, Alphabet alphabet) :
                content_(content),
                alphabet_(alphabet) {};

        STD_SEQUENCE_PROTO &operator[] (lensq index) {
            return content_[index];
        }

        const STD_SEQUENCE_PROTO &operator[] (lensq index) const {
            return content_[index];
        }

        [[nodiscard]] lensq length() const {
            return content_.size();
        }

        [[nodiscard]] AlphabetType alphabet() const {
            return alphabet_;
        }

        template<typename TYPE_OUT>
        TYPE_OUT pack() const {
            return sqapply<StdSqProto, TYPE_OUT, AlphabetType>(*this, ops::OperationPack<SequenceType, typename TYPE_OUT::SequenceType, AlphabetType>());
        }

        Rcpp::List exportToR() {
            content_.attr("alphabet") = alphabet_;
            return content_;
        }
    };
}

#endif //TIDYSQ_TYPES_STDSQPROTO_H
