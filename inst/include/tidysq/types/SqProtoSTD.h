#ifndef TIDYSQ_SQPROTOSTD_H
#define TIDYSQ_SQPROTOSTD_H

#include "SequenceProtoSTD.h"
#include "AlphabetSTD.h"
#include "../impl_Operation.h"
#include "../sqapply.h"

namespace tidysq {
    template<ProtoType PROTO>
    class SqProto<STD, PROTO> {
        std::vector<SequenceProto<STD, PROTO>> content_;
        Alphabet<STD> alphabet_;
    public:
        typedef SequenceProto<STD, PROTO> SequenceType;
        typedef Alphabet<STD> AlphabetType;

        SqProto(std::vector<SequenceProto<STD, PROTO>> content, Alphabet<STD> alphabet) :
                content_(content),
                alphabet_(alphabet) {};

        SequenceType &operator[] (lensq index) {
            return content_[index];
        }

        const SequenceType &operator[] (lensq index) const {
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
            return sqapply<SqProto<STD, PROTO>, TYPE_OUT, AlphabetType>(*this, ops::OperationPack<SequenceType, typename TYPE_OUT::SequenceType, AlphabetType>());
        }

        Rcpp::List exportToR() {
            content_.attr("alphabet") = alphabet_;
            return content_;
        }
    };
}

#endif //TIDYSQ_SQPROTOSTD_H
