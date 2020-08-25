#ifndef TIDYSQ_SQPROTOSTD_H
#define TIDYSQ_SQPROTOSTD_H

#include <utility>

#include "SequenceProtoSTD.h"
#include "Alphabet.h"
#include "../ops/OperationPack.h"
#include "../sqapply.h"

namespace tidysq {
    template<ProtoType PROTO>
    class SqProto<STD, PROTO> {
        std::vector<SequenceProto<STD, PROTO>> content_;
        Alphabet alphabet_;
    public:
        typedef SequenceProto<STD, PROTO> SequenceType;

        SqProto(const std::vector<SequenceProto<STD, PROTO>> &content, Alphabet alphabet) :
                content_(content),
                alphabet_(std::move(alphabet)) {};

        SqProto(lensq length, Alphabet alphabet) :
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

        [[nodiscard]] const Alphabet &alphabet() const {
            return alphabet_;
        }

        template<InternalType INTERNAL_OUT>
        Sq<INTERNAL_OUT> pack() const {
            return sqapply<SqProto<STD, PROTO>, Sq<INTERNAL_OUT>>(*this, ops::OperationPack<STD, PROTO, INTERNAL_OUT>());
        }
    };
}

#endif //TIDYSQ_SQPROTOSTD_H
