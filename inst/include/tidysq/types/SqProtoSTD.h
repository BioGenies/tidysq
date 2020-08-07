#ifndef TIDYSQ_SQPROTOSTD_H
#define TIDYSQ_SQPROTOSTD_H

#include <utility>

#include "SequenceProtoSTD.h"
#include "AlphabetSTD.h"
#include "../ops/OperationPack.h"
#include "../sqapply.h"

namespace tidysq {
    template<ProtoType PROTO>
    class SqProto<STD, PROTO> {
        std::vector<SequenceProto<STD, PROTO>> content_;
        Alphabet<STD> alphabet_;
    public:
        typedef SequenceProto<STD, PROTO> SequenceType;
        typedef Alphabet<STD> AlphabetType;

        SqProto(const std::vector<SequenceProto<STD, PROTO>> &content, Alphabet<STD> alphabet) :
                content_(content),
                alphabet_(std::move(alphabet)) {};

        SqProto(lensq length, Alphabet<STD> alphabet) :
                content_(std::vector<SequenceProto<STD, PROTO>>(length)),
                alphabet_(std::move(alphabet)) {};

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

        template<InternalType INTERNAL_OUT>
        Sq<INTERNAL_OUT> pack() const {
            return sqapply<SqProto<STD, PROTO>, Sq<INTERNAL_OUT>, AlphabetType>(*this, ops::OperationPack<STD, PROTO, INTERNAL_OUT>());
        }

        Rcpp::List exportToR() {
            content_.attr("alphabet") = alphabet_;
            return content_;
        }
    };
}

#endif //TIDYSQ_SQPROTOSTD_H
