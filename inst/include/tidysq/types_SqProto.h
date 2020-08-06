#ifndef TIDYSQ_TYPES_SQPROTO_H
#define TIDYSQ_TYPES_SQPROTO_H

#include "types_Sequence.h"
#include "impl_Operation.h"
#include "sqapply.h"

namespace tidysq {
    template<typename SEQUENCE_PROTO>
    class SqProto {
        std::vector<SEQUENCE_PROTO> content_;
        Alphabet alphabet_;
    public:
        typedef SEQUENCE_PROTO SequenceType;
        typedef Alphabet AlphabetType;

        SqProto(std::vector<SEQUENCE_PROTO> content, Alphabet alphabet) :
                content_(content),
                alphabet_(alphabet) {};

        SEQUENCE_PROTO &operator[] (lensq index) {
            return content_[index];
        }

        const SEQUENCE_PROTO &operator[] (lensq index) const {
            return content_[index];
        }

        lensq length() const {
            return content_.size();
        }

        AlphabetType alphabet() const {
            return alphabet_;
        }

        template<typename TYPE_OUT>
        TYPE_OUT pack() const {
            return sqapply<SqProto, TYPE_OUT, AlphabetType>(*this, ops::OperationPack<SequenceType, typename TYPE_OUT::SequenceType, AlphabetType>());
        }

        Rcpp::List exportToR() {
            content_.attr("alphabet") = alphabet_;
            return content_;
        }
    };
}

#endif //TIDYSQ_TYPES_SQPROTO_H
